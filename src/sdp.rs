//! Create and manipulate semi-definite problems.

use crate::algebra::*;
use crate::expr::Expr;
use crate::flag::Flag;
use crate::operator::*;
use crate::sdpa::*;
use crate::tools::csdp;
use crate::tools::csdp_minimize_certificate;
use arrayvec::ArrayVec;
use ndarray_linalg::{Cholesky, Eigh, UPLO};
use num::*;
use sprs::{CsMat, TriMat};
use std::ops::Index;

use ndarray::{Array1, Array2, ScalarOperand};
use std::fmt;
use std::fmt::Display;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::{AddAssign, DivAssign, Neg};
use std::path::PathBuf;

use log::*;

/// The optimization problems over flags are translated into a
/// sdp problem in the sdpa format.
///
/// Shape of the matrices:
///
/// For each i in ineqs (where i is itself a vector of inequalities):
///    A diagonal block of size i.len()
///
/// For each cs:
///    A block with the size od `cs.input_matrix`

/// An optimization problem expressed in flags algebra.
#[derive(Debug, Clone)]
pub struct Problem<N, F: Flag> {
    /// Set of contraint inequalities.
    pub ineqs: Vec<Ineq<N, F>>,
    /// Set of Cauchy-Schwarz inequalities to be used.
    pub cs: Vec<MulAndUnlabel<F>>,
    /// Vector to be optimized.
    pub obj: QFlag<N, F>,
}

// Write a line of a sdpa file
fn write_coeff<N>(
    file: &mut BufWriter<File>,
    mat_num: usize,
    block_num: usize,
    i: usize,
    j: usize,
    value: N,
) -> io::Result<()>
where
    N: Display + Zero + PartialEq,
{
    if value != N::zero() {
        writeln!(
            file,
            "{} {} {} {} {}",
            mat_num,
            block_num + 1,
            i + 1,
            j + 1,
            value
        )?;
    }
    Ok(())
}

// Write a coefficient `value` on position (i ,i) for an inequality,
// or `value` and `-value` on positions (2i, 2i), (2i+1, 2i+1) for an equality
fn write_ineq_coeff<N>(
    file: &mut BufWriter<File>,
    mat_num: usize,
    block_num: usize,
    i: usize,
    value: N,
    equality: bool,
) -> io::Result<()>
where
    N: Display + Zero + PartialEq + Neg<Output = N> + Copy,
{
    if equality {
        write_coeff(file, mat_num, block_num, 2 * i, 2 * i, value)?;
        write_coeff(file, mat_num, block_num, 2 * i + 1, 2 * i + 1, -value)?;
    } else {
        write_coeff(file, mat_num, block_num, i, i, value)?;
    }
    Ok(())
}

impl<N, F: Flag> Problem<N, F> {
    /// Create a minimization problem with the argument as objective function.
    pub fn minimize(obj: QFlag<N, F>) -> Self
    where
        N: Display + Num + FromPrimitive + Clone + Neg<Output = N>,
    {
        Self {
            ineqs: vec![
                flags_are_nonnegative(obj.basis),
                total_sum_is_one(obj.basis),
            ],
            cs: obj.basis.all_cs(),
            obj,
        }
    }
    /// Panic if the size of the basis involved are inconsistent.
    pub fn check(&self) {
        let b = self.obj.basis;
        for ineq in &self.ineqs {
            assert_eq!(ineq.meta.basis, b);
        }
        for cs in &self.cs {
            assert_eq!(cs.output_basis(), b);
        }
    }
    pub(crate) fn view<'a>(&'a self, selector: &'a Selector) -> ProblemView<'a, N, F> {
        self.check();
        let ineqs = Select {
            selected: &self.ineqs,
            selector: &selector.ineqs,
        };
        let cs = Select {
            selected: &self.cs,
            selector: &selector.cs,
        };
        ProblemView {
            obj: &self.obj,
            ineqs,
            cs,
            cs_subspace: selector.cs_subspace_constraints().collect(),
        }
    }
}

fn format_count(count: usize, max: usize) -> String {
    if count == max {
        format!("{count}")
    } else {
        format!("{count}/{max}")
    }
}

impl<'a, N, F> ProblemView<'a, N, F>
where
    N: Display + Neg<Output = N> + Zero + Copy + PartialEq,
    F: Flag,
{
    /// Write a description of the problem as a comment in the .sdpa file
    fn write_header<W: Write>(&self, file: &mut W) -> io::Result<()> {
        writeln!(file, "* Semi-Definite Problem generated by Rust")?;
        writeln!(file, "*")?;
        writeln!(file, "* Flags: {}", F::NAME)?;
        writeln!(
            file,
            "* Basis: {:} ({} flags)",
            self.obj.basis,
            self.obj.basis.get().len()
        )?;
        writeln!(
            file,
            "* {} groups of linear constraints ({})",
            self.ineqs.len(),
            format_count(
                self.ineqs.iter().map(|i| { i.len() }).sum::<usize>(),
                self.ineqs
                    .iter()
                    .map(|ineq| { ineq.selected.data.len() })
                    .sum::<usize>()
            )
        )?;
        writeln!(
            file,
            "* {} Cauchy-Schwarz constraints",
            format_count(self.cs.len(), self.cs.selected.len())
        )?;
        writeln!(
            file,
            "* {} additional constraints on Cauchy-Schwarz matrices",
            self.cs_subspace.len(),
        )?;
        writeln!(file, "*")?;
        write!(file, "* Minimizing: {}", self.obj.expr)?;
        match self.obj.scale {
            1 => writeln!(file)?,
            s => writeln!(file, " (scale {s})")?,
        }
        writeln!(file, "* Under the constraints:")?;
        for ineqs in self.ineqs.iter() {
            writeln!(
                file,
                "* # {} ({})",
                ineqs.meta(),
                format_count(ineqs.len(), ineqs.selected.data.len())
            )?;
        }
        for cs in self.cs.iter() {
            writeln!(
                file,
                "* # Cauchy-Schwarz: {}x{}; {:?}",
                cs.1.split.left_basis().size,
                cs.1.split.right_basis().size,
                cs.1.split.right_basis().t,
            )?;
        }
        writeln!(file, "*")?;
        Ok(())
    }
    /// Write the semi-definite program in the file `filename` in the sdpa format.
    pub fn write_sdpa(&self, filename: &str) -> io::Result<()> {
        let mut filename = PathBuf::from(filename);
        let _ = filename.set_extension("sdpa");
        info!("Writing problem in {}", filename.display());
        let mut file = BufWriter::new(File::create(&filename)?);
        self.write_header(&mut file)?;
        debug!("Generating Cauchy-Schwarz inequalities");
        let cs_mat: Vec<_> = self.cs.get();
        writeln!(file, "{}", self.obj.data.len() + self.cs_subspace.len())?;
        // Line 2: Number of blocks (one for each constraint)
        writeln!(file, "{}", self.ineqs.len() + self.cs.len())?;
        // Line 3: Sizes of the blocks
        for ineq in self.ineqs.iter() {
            assert!(ineq.len() > 0);
            write!(file, "-{} ", ineq.len_spliting_equalities())?;
        }
        for split in &cs_mat {
            write!(file, "{} ", split[0].rows())?;
        }
        writeln!(file)?;
        // Line 4: vector ai
        // ai is the needed coefficient for the flag Fi
        for v in &self.obj.data {
            write!(file, "{v} ")?;
        }
        for _ in &self.cs_subspace {
            write!(file, "0 ")?;
        }
        write!(file, "\n\n")?;
        // Lines 5+: body
        // Matrix 0: Objective
        for (block_num, ineq) in self.ineqs.iter().enumerate() {
            for (i, ineq_data) in ineq.iter().enumerate() {
                write_ineq_coeff(
                    &mut file,
                    0,
                    block_num,
                    i,
                    ineq_data.bound,
                    ineq.meta().equality,
                )?;
            }
        }
        writeln!(file)?;
        // Matrices 1+:
        // Inequaltity blocks
        for (block_num, ineq) in self.ineqs.iter().enumerate() {
            for (i, ineq_data) in ineq.iter().enumerate() {
                for (mat_num, &v) in ineq_data.flag.iter() {
                    write_ineq_coeff(
                        &mut file,
                        mat_num + 1,
                        block_num,
                        i,
                        v,
                        ineq.meta().equality,
                    )?;
                }
            }
        }
        writeln!(file)?;
        // Cs blocks
        let block_offset = self.ineqs.len();
        for (block_num, line) in cs_mat.iter().enumerate() {
            for (mat_num, matrix) in line.iter().enumerate() {
                for (&v, (i, j)) in matrix {
                    if i <= j {
                        write_coeff(&mut file, mat_num + 1, block_num + block_offset, i, j, v)?;
                    }
                }
            }
        }
        // Cs subspace constraints
        let mat_offset = self.obj.data.len();
        for (i_mat, (block_num, matrix)) in self.cs_subspace.iter().enumerate() {
            for (&v, (i, j)) in *matrix {
                if i <= j {
                    let mat_num = mat_offset + i_mat + 1;
                    write_coeff(&mut file, mat_num, block_num + block_offset, i, j, v)?;
                }
            }
        }
        Ok(())
    }
}

impl<N, F> Problem<N, F>
where
    N: Display + Zero + Copy + PartialEq + Neg<Output = N>,
    F: Flag,
{
    /// Write the semi-definite program in the file `filename` in the sdpa format.
    pub fn write_sdpa(&self, filename: &str) -> io::Result<()> {
        self.view(&Selector::new(self)).write_sdpa(filename)
    }
    /// Rescale the objective according to its scale field.
    /// If this method is not used, the output of the sdp solver may need to be rescaled.
    pub fn no_scale(mut self) -> Self
    where
        N: DivAssign + ScalarOperand + FromPrimitive,
    {
        self.obj = self.obj.no_scale();
        self
    }
    /// Solve the sdp using the CSDP solver.
    pub fn solve_csdp(&self, filename: &str) -> Result<f64, Error> {
        self.write_sdpa(filename)?;
        self.run_csdp(filename, None, false)
    }
    pub fn run_csdp(
        &self,
        name: &str,
        initial_solution: Option<&str>,
        minimize_certificate: bool,
    ) -> Result<f64, Error> {
        let filename = format!("{name}.sdpa");
        match if minimize_certificate {
            csdp_minimize_certificate(&filename, initial_solution)
        } else {
            csdp(&filename, initial_solution)
        } {
            Ok(v) => Ok(v / self.obj.scale as f64),
            e => e,
        }
    }
}

// Selectors

/// An object that specify a subset of constraint of a flag problem.
#[derive(Debug, Clone, PartialEq)]
pub struct Selector {
    pub ineqs: Vec<Vec<usize>>,
    cs: Vec<VecCsMode>,
    cs_subspace: Vec<Subspace>,
}

#[derive(Debug, Clone, PartialEq)]
struct Subspace {
    basis_len: usize,
    classes: usize,
    orth: Vec<(CSMode, CsMat<f64>)>,
}

impl Subspace {
    fn new<F: Flag>(cs: &MulAndUnlabel<F>) -> Self {
        Self {
            basis_len: Unlabel::total(cs.split.left_basis()).get().0.len(),
            classes: cs.invariant_classes().get().into_iter().max().unwrap() + 1,
            orth: Vec::new(),
        }
    }
    fn dim(&self, mode: CSMode) -> usize {
        match mode {
            Simple => self.basis_len,
            Invariant => self.classes,
            AntiInvariant => self.basis_len - self.classes,
        }
    }
    fn restrict(&mut self, mode: CSMode, mat: CsMat<f64>) {
        self.orth.push((mode, mat))
    }
    fn orthogonal_matrices(&self, mode: CSMode) -> impl Iterator<Item = &CsMat<f64>> {
        self.orth
            .iter()
            .filter(move |(mode2, _)| *mode2 == mode)
            .map(|(_, mat)| mat)
    }
}

/// Specifies a symmetry reduction for a product-and-unlabel matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSMode {
    Simple,
    Invariant,
    AntiInvariant,
}

use CSMode::*;

type VecCsMode = ArrayVec<CSMode, 2>;

/// Identifies a product-and-unlabel matrix with a symmetry reduction
#[derive(Debug, Clone)]
pub struct CauchySchwarzMatrix<F: Flag>(pub CSMode, pub MulAndUnlabel<F>);

impl<F: Flag> Display for CauchySchwarzMatrix<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.0 {
            Simple => write!(f, "{}", self.1),
            Invariant => write!(f, "{} invariant", self.1),
            AntiInvariant => write!(f, "{} anti-invariant", self.1),
        }
    }
}

impl<F: Flag> CauchySchwarzMatrix<F> {
    pub fn get(&self) -> Vec<CsMat<i64>> {
        match self.0 {
            Simple => self.1.get(),
            Invariant => self.1.reduced().get().0,
            AntiInvariant => self.1.reduced().get().1,
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Select<'a, A, B> {
    selected: &'a A,
    selector: &'a B,
}

#[derive(Debug, Clone)]
pub(crate) struct SelectIter<'a, A, B> {
    content: &'a Select<'a, A, B>,
    iter: std::ops::Range<usize>,
}

/// A reference to a `Problem` filtered by a `Selector`.
#[derive(Debug, Clone)]
pub struct ProblemView<'a, N, F: Flag> {
    pub obj: &'a QFlag<N, F>,
    pub(crate) ineqs: IneqsSelect<'a, N, F>,
    pub(crate) cs: CsSelect<'a, F>,
    pub(crate) cs_subspace: Vec<(usize, &'a CsMat<f64>)>,
}

// Sub-selector types
type CsSelect<'a, F> = Select<'a, Vec<MulAndUnlabel<F>>, Vec<VecCsMode>>;
type IneqsSelect<'a, N, F> = Select<'a, Vec<Ineq<N, F>>, Vec<Vec<usize>>>;
type IneqSelect<'a, N, F> = Select<'a, Ineq<N, F>, Vec<usize>>;

#[derive(Debug, Clone)]
pub struct CsIter<'a, F: Flag> {
    content: CsSelect<'a, F>,
    iter: std::ops::Range<usize>,
    next: Option<CauchySchwarzMatrix<F>>,
}

type IneqsIter<'a, N, F> = SelectIter<'a, Vec<Ineq<N, F>>, Vec<Vec<usize>>>;
type IneqIter<'a, N, F> = SelectIter<'a, Ineq<N, F>, Vec<usize>>;

impl<'a, F: Flag> CsSelect<'a, F> {
    pub fn iter(&self) -> CsIter<'a, F> {
        CsIter {
            content: self.clone(),
            iter: 0..self.selector.len(),
            next: None,
        }
    }
    pub fn len(&self) -> usize {
        self.iter().count()
    }
    pub fn get(&self) -> Vec<Vec<CsMat<i64>>>
    where
        F: Flag,
    {
        self.iter().map(|mat| mat.get()).collect()
    }
}

impl<'a, N, F: Flag> Index<usize> for IneqSelect<'a, N, F> {
    type Output = IneqData<N>;

    fn index(&self, i: usize) -> &Self::Output {
        &self.selected.data[self.selector[i]]
    }
}

impl<'a, N, F: Flag> IneqsSelect<'a, N, F> {
    pub fn iter(&'a self) -> IneqsIter<'a, N, F> {
        SelectIter {
            content: self,
            iter: 0..self.selector.len(),
        }
    }
    /// Compute the number of group of inequalities
    pub fn len(&self) -> usize {
        self.selector
            .iter()
            .filter(|select| !select.is_empty())
            .count()
    }
    pub fn get(&self, i: usize) -> Option<IneqSelect<'a, N, F>> {
        if self.selector[i].is_empty() {
            None
        } else {
            Some(Select {
                selector: &self.selector[i],
                selected: &self.selected[i],
            })
        }
    }
}

impl<'a, N, F: Flag> IneqSelect<'a, N, F> {
    pub fn iter(&'a self) -> IneqIter<'a, N, F> {
        SelectIter {
            content: self,
            iter: 0..self.selector.len(),
        }
    }
    /// Number of (in)equalities in the selected group.
    pub fn len(&self) -> usize {
        self.selector.len()
    }
    /// Number of equalities in the selected group where equalities count for 2 (for ≥ and ≤).
    pub fn len_spliting_equalities(&self) -> usize {
        let len = self.len();
        if self.meta().equality {
            len * 2
        } else {
            len
        }
    }
    pub fn meta(&self) -> &IneqMeta<N, F> {
        &self.selected.meta
    }
}

// Iterators
impl<'a, F: Flag> Iterator for CsIter<'a, F> {
    type Item = CauchySchwarzMatrix<F>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(x) = self.next.take() {
            Some(x)
        } else {
            // FIXME
            match self.iter.next() {
                None => None,
                Some(i) => {
                    let cs = self.content.selected[i];
                    let mat = |m| CauchySchwarzMatrix(m, cs);
                    let select_i = &self.content.selector[i];
                    match select_i.len() {
                        0 => self.next(),
                        1 => Some(mat(select_i[0])),
                        2 => {
                            self.next = Some(mat(select_i[0]));
                            Some(mat(select_i[1]))
                        }
                        _ => unimplemented!(),
                    }
                }
            }
        }
    }
}

impl<'a, N, F: Flag> Iterator for IneqIter<'a, N, F> {
    type Item = &'a IneqData<N>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|i| &self.content[i])
    }
}

impl<'a, N, F: Flag> Iterator for IneqsIter<'a, N, F> {
    type Item = IneqSelect<'a, N, F>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            Some(i) => match self.content.get(i) {
                None => self.next(),
                some => some,
            },
            None => None,
        }
    }
}

type Id = (usize, CSMode);

impl Selector {
    pub fn new<N, F: Flag>(prob: &Problem<N, F>) -> Self {
        let mut simple = ArrayVec::new();
        simple.push(Simple); //FIXME
                             //simple.push(Invariant);
        (Self {
            ineqs: prob
                .ineqs
                .iter()
                .map(|ineq| (0..ineq.data.len()).collect())
                .collect(),
            cs: vec![simple; prob.cs.len()],
            cs_subspace: prob.cs.iter().map(Subspace::new).collect(), // Can be expensive
        })
        .variant_reduced(prob) // FIXME
    }
    pub fn variant_reduced<N, F>(mut self, problem: &Problem<N, F>) -> Self
    where
        F: Flag,
    {
        for (mode, cs) in self.cs.iter_mut().zip(&problem.cs) {
            if mode.as_slice().contains(&Simple) {
                mode.clear();
                mode.push(Invariant);
                let class = cs.invariant_classes().get();
                let n = class.len();
                if class[n - 1] < n - 1 {
                    mode.push(AntiInvariant);
                }
            }
        }
        self
    }
    pub fn cs_subspace_constraints(&self) -> impl Iterator<Item = (usize, &CsMat<f64>)> {
        self.cs
            .iter()
            .zip(self.cs_subspace.iter())
            .flat_map(|(cs_modes, subspace)| {
                cs_modes
                    .iter()
                    .map(move |&mode| subspace.orthogonal_matrices(mode))
            })
            .enumerate()
            .flat_map(|(i, iter_on_cs)| iter_on_cs.map(move |mat| (i, mat)))
    }
    pub fn weight(&self) -> (usize, usize) {
        let count_ineq = self.ineqs.iter().map(|l| l.len()).sum();
        (count_ineq, self.cs.len())
    }
    pub fn refine_with_certificate(&self, cert: &Certificate<f64>, protected: &[bool]) -> Self {
        assert_eq!(self.ineqs.len(), protected.len());
        let mut res = self.clone();
        for (i, protect) in protected.iter().enumerate() {
            if !protect {
                for j in (0..self.ineqs[i].len()).rev() {
                    if match cert.x[i].get(j, j) {
                        None => true,
                        Some(v) => v < &1e-12f64,
                    } {
                        let _ = res.ineqs[i].remove(j);
                    }
                }
            }
        }
        res
    }
    pub fn cs_dim(&self, (i, mode): Id) -> usize {
        self.cs_subspace[i].dim(mode)
    }
    pub fn remove_cs(&self, (i, mode): Id) -> Result<Self, String> {
        match self.cs[i].iter().position(|&m| m == mode) {
            None => Err(format!("No mode {mode:?} for index {i} to remove")),
            Some(j) => {
                let mut res = self.clone();
                let _ = res.cs[i].swap_remove(j);
                Ok(res)
            }
        }
    }
    pub fn cs_vec(&self) -> Vec<Id> {
        let mut res = Vec::new();
        for (i, vec) in self.cs.iter().enumerate() {
            for &mode in vec.as_slice() {
                res.push((i, mode))
            }
        }
        res
    }
    pub fn restrict_cs(&self, (i, mode): Id, mat: CsMat<f64>) -> Result<Self, String> {
        if self.cs[i].iter().any(|&m| m == mode) {
            // If the id is valid
            let mut res = self.clone();
            res.cs_subspace[i].restrict(mode, mat);
            Ok(res)
        } else {
            Err(format!("No mode {mode:?} for index {i} found"))
        }
    }
}

/// A certificate to a sdp problem as given by CSDP
#[derive(Debug, Clone)]
pub struct Certificate<N> {
    pub y: Vec<N>,
    pub z: Vec<CsMat<N>>, // matrix 1 in sdpa file
    pub x: Vec<CsMat<N>>, // matrix 2 in sdpa file
}

impl Certificate<f64> {
    pub fn from_file_select<N, F: Flag>(
        pb: &ProblemView<'_, N, F>,
        name: &str,
    ) -> io::Result<Self> {
        let file = File::open(name)?;
        let mut buf = BufReader::new(file).lines();
        // The value of the vector y is on the first line
        let y = buf
            .next()
            .unwrap()?
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect();
        // Prepare the space for the matrices z and x
        let mut tri_mat = [Vec::new(), Vec::new()];
        for ineq in pb.ineqs.iter() {
            let len = ineq.len_spliting_equalities(); // Bigger matrices to gather the coefficients
            assert!(len > 0);
            for m in &mut tri_mat {
                m.push(TriMat::new((len, len)))
            }
        }
        for cs in &pb.cs.get() {
            // We can avoid some memory usage here
            let n = cs[0].rows();
            for m in &mut tri_mat {
                m.push(TriMat::new((n, n)))
            }
        }
        for line in buf {
            let l = line.unwrap().parse::<SdpaCoeff>().unwrap();
            tri_mat[l.mat - 1][l.block - 1].add_triplet(l.i - 1, l.j - 1, l.val);
            // If the value is not on the diagonal, add the symmetric coeff
            if l.i != l.j {
                tri_mat[l.mat - 1][l.block - 1].add_triplet(l.j - 1, l.i - 1, l.val);
            }
        }
        // condense inequality matrices
        for (i, ineq) in pb.ineqs.iter().enumerate() {
            if ineq.meta().equality {
                for blocks in &mut tri_mat {
                    let old_mat = &mut blocks[i];
                    let len = ineq.len();
                    let mut new_mat = TriMat::with_capacity((len, len), old_mat.nnz());
                    for (&val, (i, j)) in old_mat.triplet_iter() {
                        assert_eq!(i, j);
                        let new_val = if i % 2 == 0 { val } else { -val };
                        new_mat.add_triplet(i / 2, i / 2, new_val)
                    }
                    *old_mat = new_mat
                }
            }
        }
        // convert the triplet matrices to sparse matrices
        let mut sprs_mat = tri_mat
            .iter()
            .map(|mat| mat.iter().map(|block| block.to_csc()).collect());
        Ok(Self {
            y,
            z: sprs_mat.next().unwrap(), // First matrix: z
            x: sprs_mat.next().unwrap(), // Second matrix: x
        })
    }
    pub fn to_file(&self, name: &str) -> io::Result<()> {
        let mut w = BufWriter::new(File::create(name)?);
        // The value of the vector y is on the first line
        for v in &self.y {
            write!(w, "{v} ")?;
        }
        writeln!(w)?;
        for (num_mat, matrix) in [&self.z, &self.x].iter().enumerate() {
            for (block, mat) in matrix.iter().enumerate() {
                for (v, (i, j)) in mat {
                    if i <= j {
                        writeln!(w, "{} {} {} {} {}", num_mat + 1, block + 1, i + 1, j + 1, v)?;
                    }
                }
            }
        }
        Ok(())
    }
    pub fn value_primal<F, N>(&self, pb: &Problem<N, F>) -> f64
    where
        N: Clone + ToPrimitive,
        F: Flag,
    {
        let mut res = 0.;
        assert_eq!(pb.obj.data.len(), self.y.len());
        for (ai, yi) in pb.obj.data.iter().zip(self.y.iter()) {
            let ai: f64 = NumCast::from(ai.clone()).unwrap();
            res += yi * ai;
        }
        res
    }
    // C * X, here C is 0 for Cauchy-Schwarz inequalities
    pub fn value_dual<F, N>(&self, pb: &Problem<N, F>) -> f64
    where
        N: Clone + ToPrimitive,
        F: Flag,
    {
        let mut res = 0.;
        // For each block of inequalities
        for (block, ineqs) in pb.ineqs.iter().enumerate() {
            for (v, (i, j)) in &self.x[block] {
                assert_eq!(i, j);
                let bound: f64 = NumCast::from(ineqs.data[i].bound.clone()).unwrap();
                res += v * bound;
            }
        }
        res
    }
    pub fn values<F, N>(&self, pb: &Problem<N, F>) -> (f64, f64)
    where
        N: Clone + ToPrimitive,
        F: Flag,
    {
        (self.value_primal(pb), self.value_dual(pb))
    }
    pub fn to_vec<F>(&self, b: Basis<F>, threshold: f64) -> QFlag<f64, F>
    where
        F: Flag,
    {
        b.qflag_from_vec(
            self.y
                .iter()
                .map(|x| if x.abs() < threshold { 0. } else { *x })
                .collect(),
        )
    }
    pub fn with_threshold(mut self, threshold: f64) -> Self {
        for x in &mut self.y {
            if x.abs() < threshold {
                *x = 0.
            }
        }
        self
    }
    pub fn diag_coeffs(&self, block: usize, n: usize) -> Vec<f64> {
        assert_eq!(self.x[block].cols(), n);
        assert_eq!(self.x[block].rows(), n);
        let mut res = vec![0.; n];
        for (&v, (i, j)) in &self.x[block] {
            assert_eq!(i, j);
            res[i] += v;
        }
        res
    }
}

pub(crate) fn condense<N, F>(ineq: IneqSelect<'_, N, F>, coeff: &[N]) -> (QFlag<N, F>, N)
where
    N: Num + Clone + ScalarOperand + AddAssign,
    F: Flag,
{
    assert_eq!(ineq.len(), coeff.len());
    assert!(!coeff.is_empty());
    let mut bound = N::zero();
    let mut res: Array1<N> = Array1::zeros(ineq.selected.data[0].flag.dim());
    for (c, ineq_data) in coeff.iter().zip(ineq.iter()) {
        if c != &N::zero() {
            for (i, val) in ineq_data.flag.iter() {
                res[i] += val.clone() * c.clone();
            }
            bound += ineq_data.bound.clone() * c.clone()
        }
    }
    (
        QFlag {
            basis: ineq.meta().basis,
            data: res,
            scale: 1,
            expr: ineq.meta().flag_expr.clone(),
        },
        bound,
    )
}

fn hadamard<N>(dense: &Array2<N>, sprs: &CsMat<i64>) -> N
where
    N: Num + Clone + FromPrimitive,
{
    let mut res = N::zero();
    for (&v, (i, j)) in sprs {
        res = res + N::from_i64(v).unwrap() * dense[(i, j)].clone()
    }
    res
}

pub fn condense_cs<N, F>(cs: &CauchySchwarzMatrix<F>, cert: &Array2<N>) -> QFlag<N, F>
where
    N: Num + Clone + FromPrimitive,
    F: Flag,
{
    let mut data = Vec::new();
    for m in &cs.get() {
        data.push(hadamard(cert, m));
    }
    QFlag {
        basis: cs.1.output_basis(),
        data: Array1::from(data),
        scale: 1,
        expr: Expr::unknown("Cauchy-Schwarz".into()),
    }
}

pub fn cholesky_qflags<F>(cs: &CauchySchwarzMatrix<F>, mat: &Array2<f64>) -> Vec<QFlag<f64, F>>
where
    F: Flag,
{
    let mut res = Vec::with_capacity(mat.ncols());

    let cholesky0 = Cholesky::cholesky(mat, UPLO::Upper).unwrap();
    // Fixme

    let cholesky = match cs.0 {
        Invariant => {
            let inv_mat = crate::density::class_matrices(&cs.1.invariant_classes().get())
                .0
                .map(|&x| x as f64);
            (&inv_mat * &cholesky0.t()).t().to_owned()
        }
        AntiInvariant => {
            let inv_mat = crate::density::class_matrices(&cs.1.invariant_classes().get())
                .1
                .map(|&x| x as f64);
            (&inv_mat * &cholesky0.t()).t().to_owned()
        }
        _ => cholesky0,
    };

    //FIXME
    for i in 0..cholesky.nrows() {
        let data = cholesky.row(i).to_owned();
        assert_eq!(data.len(), cs.1.split.left_basis().get().len());
        res.push(QFlag {
            basis: cs.1.split.left_basis(),
            data,
            scale: 1,
            expr: Expr::unknown("Cholesky".into()),
        })
    }

    res
}

pub fn eigenvectors_qflags<F>(
    cs: &CauchySchwarzMatrix<F>,
    mat: &Array2<f64>,
) -> Vec<(f64, QFlag<f64, F>)>
where
    F: Flag,
{
    // FIXME: code duplication
    let mut res = Vec::with_capacity(mat.ncols());

    let (eigenvalues, eigenvectors) = mat.eigh(UPLO::Lower).unwrap();

    let eigenvectors = eigenvectors.t().to_owned();

    let eigenvectors = match cs.0 {
        Invariant => {
            let inv_mat = crate::density::class_matrices(&cs.1.invariant_classes().get())
                .0
                .map(|&x| x as f64);
            (&inv_mat * &eigenvectors.t()).t().to_owned()
        }
        AntiInvariant => {
            let inv_mat = crate::density::class_matrices(&cs.1.invariant_classes().get())
                .1
                .map(|&x| x as f64);
            (&inv_mat * &eigenvectors.t()).t().to_owned()
        }
        Simple => eigenvectors,
    };

    for i in 0..eigenvectors.nrows() {
        let data = eigenvectors.row(i).to_owned();
        assert_eq!(data.len(), cs.1.split.left_basis().get().len());
        res.push((
            eigenvalues[i],
            QFlag {
                basis: cs.1.split.left_basis(),
                data,
                scale: 1,
                expr: Expr::unknown("Eigenvectors".into()),
            },
        ))
    }

    res
}
