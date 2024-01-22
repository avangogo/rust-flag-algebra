//! Computing and stocking operators of the flag algebra.

use crate::algebra::QFlag;
use crate::combinatorics::*;
use crate::density::*;
use crate::expr::Expr;
use crate::flag::Flag;
use log::*;
use ndarray::*;
use num::*;
use serde::de::DeserializeOwned;
use serde::Serialize;
use sprs::CsMat;
use std::error::Error;
use std::fmt;
use std::fmt::{Display, Formatter};
use std::fs;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::marker::PhantomData;
use std::ops::*;
use std::path::*;
use std::rc::Rc;

/// A trait for flag operators that
/// can be saved in a file once computed
/// for the first time.
///
/// `A` is the type of the stored object.
/// The operator operate on flags of type `F`.
pub trait Savable<A, F>
where
    A: Serialize + DeserializeOwned,
    F: Flag,
{
    /// Name of the file where the operator can be saved.
    fn filename(&self) -> String;
    /// Compute the object.
    fn create(&self) -> A;
    /// Path to the corresponding file.
    fn file_path(&self) -> PathBuf {
        let mut filename = PathBuf::from("./data");
        filename.push(Path::new(F::NAME));
        filename.push(self.filename());
        let _ = filename.set_extension("dat");
        filename
    }
    /// (Re)create the object, save it in the corresponding file and return it.
    fn create_and_save(&self, path: &Path) -> A {
        info!("Creating {}", path.display());
        let value = self.create();
        let file = File::create(path).unwrap();
        let buf = BufWriter::new(file);
        bincode::serialize_into(buf, &value).unwrap();
        value
    }
    /// Load the object if the file exists and is valid.
    fn load(&self, path: &Path) -> Result<A, Box<dyn Error>> {
        let file = File::open(path)?;
        let mut buf = BufReader::new(file);
        let data = bincode::deserialize_from(&mut buf)?;
        Ok(data)
    }
    /// Function to automatically load the object if the file exists and
    /// is valid, or create and save it otherwise.
    fn get(&self) -> A {
        let path = self.file_path();
        if path.exists() {
            debug!("Loading {}", path.display());
            match self.load(&path) {
                Ok(v) => {
                    trace!("Done");
                    v
                }
                Err(e) => {
                    error!("Failed to load {}: {}", path.display(), e);
                    self.create_and_save(&path)
                }
            }
        } else {
            let dir = path.parent().unwrap();
            match fs::create_dir_all(dir) {
                Ok(()) => self.create_and_save(&path),
                Err(e) => {
                    eprintln!("Cannot create {}.", dir.display());
                    panic!("{}", e);
                }
            }
        }
    }
}

/// Type (or root) of a flag.
/// It is identified by its size and its id in the list of flags of that size.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Clone)]
pub struct Type<F: Flag> {
    /// Size of the type.
    pub size: usize,
    /// Index of the type in the list of unlabeled flags of this size.
    pub id: usize,
    /// Retains the kind of flag (graph, ...)
    phantom: PhantomData<F>,
}

impl<F: Flag> Type<F> {
    /// Constructor for the type.
    pub fn new(size: usize, id: usize) -> Self {
        Self {
            size,
            id,
            phantom: PhantomData,
        }
    }
    /// Create a type of size 0.
    pub fn empty() -> Self {
        Self::new(0, 0)
    }
    /// Return wether the input has size zero.
    pub fn is_empty(self) -> bool {
        self == Self::empty()
    }
    /// Write a string that identifies the type.
    fn to_string_suffix(self) -> String {
        if self.is_empty() {
            String::new()
        } else {
            format!("_type_{}_id_{}", self.size, self.id)
        }
    }
    /// Create the type corresponding to g
    pub fn from_flag<G>(g: &G) -> Self
    where
        F: Flag,
        G: Into<F> + Clone,
    {
        let f: F = g.clone().into();
        let size = f.size();
        let reduced_f = f.canonical();
        let id = Basis::new(size)
            .get()
            .binary_search(&reduced_f)
            .expect("Flag not found");
        Self::new(size, id)
    }
    /// Iterate on all types of a given size.
    pub fn types_with_size(size: usize) -> impl Iterator<Item = Self>
    where
        F: Flag,
    {
        let n_types = Basis::<F>::new(size).get().len();
        (0..n_types).map(move |id| Self::new(size, id))
    }
    /// Print the type identifier in a short way.
    pub fn print_concise(self) -> String {
        if self.is_empty() {
            String::new()
        } else {
            format!("{},id{}", self.size, self.id)
        }
    }
}

//============ Basis ===========

/// Identifier for the set of flags with given size and type
/// (in the sense of a labeled subgraph).
///
/// The kind of flag is determined by the associated Rust datatype.
/// For instance `Basis<Graph>` is the Rust type for a basis of graphs.
///
/// `basis.get()` returns an ordered vector containing all corresponding flags.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Clone)]
pub struct Basis<F: Flag> {
    /// Number of vertices in the flags of the basis.
    pub size: usize,
    /// Type of the flags of the basis.
    pub t: Type<F>,
}

/// # Defining a Basis
impl<F: Flag> Basis<F> {
    /// Constructor for a basis.
    fn make(size: usize, t: Type<F>) -> Self {
        assert!(t.size <= size);
        Self { size, t }
    }
    /// Basis of flag with `size` vertices and without type.
    ///```
    /// use flag_algebra::*;
    /// use flag_algebra::flags::Graph;
    ///
    /// // Set of graphs of size 3
    /// // (the kind of flag -Graph- is deduced by type inference)
    /// let basis = Basis::new(3);
    /// let size_3_graphs: Vec<Graph> = basis.get();
    /// assert_eq!(size_3_graphs.len(), 4);
    ///
    /// // With explicit type annotation
    /// let same_basis: Basis<Graph> = Basis::new(3);
    ///```
    pub fn new(size: usize) -> Self {
        Self::make(size, Type::empty())
    }
    /// Basis of flag with same size as `self` and type `t`.
    ///```
    /// use flag_algebra::*;
    /// use flag_algebra::flags::Graph;
    ///
    /// // Basis of graphs of size 4 rooted on an edge
    /// let edge = Graph::new(2, &[(0, 1)]);
    /// let t = Type::from_flag(&edge);
    /// let basis: Basis<Graph> = Basis::new(4).with_type(t);
    ///```
    pub fn with_type(self, t: Type<F>) -> Self {
        Self::make(self.size, t)
    }
    /// Basis of flag with same size as `self` without type `t`.
    pub fn without_type(self) -> Self {
        self.with_type(Type::empty())
    }
    /// Basis of flag with `size` vertices and same type as `self`.
    pub fn with_size(self, size: usize) -> Self {
        Self::make(size, self.t)
    }
    /// Print the basis information in a short way.
    pub fn print_concise(self) -> String {
        if self.t == Type::empty() {
            format!("{}", self.size)
        } else {
            format!("{},{},id{}", self.size, self.t.size, self.t.id)
        }
    }
}

impl<F: Flag> Display for Type<F> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        if *self == Type::empty() {
            write!(f, "Empty type")
        } else {
            write!(f, "Type of size {} (id {})", self.size, self.id)
        }
    }
}

impl<F: Flag> Display for Basis<F> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        if self.t == Type::empty() {
            write!(f, "Flags of size {} without type", self.size)
        } else {
            write!(f, "Flags of size {} with type {})", self.size, self.t)
        }
    }
}

impl<F: Flag> Mul for Basis<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.t, rhs.t);
        self.with_size(self.size + rhs.size - self.t.size)
    }
}

impl<F: Flag> Div for Basis<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        assert_eq!(self.t, rhs.t);
        assert!(self.size >= rhs.size);
        self.with_size(self.size - rhs.size + self.t.size)
    }
}

impl<F: Flag> Savable<Vec<F>, F> for Basis<F> {
    fn filename(&self) -> String {
        format!("flags_{}{}", self.size, self.t.to_string_suffix())
    }
    fn create(&self) -> Vec<F> {
        if self.t == Type::empty() {
            if self.size == 0 {
                F::size_zero_flags()
            } else {
                F::generate_next(&self.with_size(self.size - 1).get())
            }
        } else {
            let type_basis = Self::new(self.t.size).get();
            let type_flag = &type_basis[self.t.id];
            F::generate_typed_up(type_flag, &self.without_type().get())
        }
    }
}

//================ Split count
/// Operator used for flag multiplication
/// .get() returns a vector of matrices M
/// where M[i][j, k] is the number of ways to split
/// i into j and k
#[derive(Debug, Clone)]
pub struct SplitCount<F: Flag> {
    left_size: usize,
    right_size: usize,
    type_: Type<F>,
}

impl<F: Flag> SplitCount<F> {
    pub fn make(left_size: usize, right_size: usize, type_: Type<F>) -> Self {
        assert!(type_.size <= left_size);
        assert!(type_.size <= right_size);
        Self {
            left_size,
            right_size,
            type_,
        }
    }

    pub fn from_input(left: &Basis<F>, right: &Basis<F>) -> Self {
        assert_eq!(left.t, right.t);
        Self::make(left.size, right.size, left.t)
    }

    pub fn left_basis(&self) -> Basis<F> {
        Basis::make(self.left_size, self.type_)
    }

    pub fn right_basis(&self) -> Basis<F> {
        Basis::make(self.right_size, self.type_)
    }

    fn output_basis(&self) -> Basis<F> {
        Basis::make(
            self.right_size + self.left_size - self.type_.size,
            self.type_,
        )
    }
    pub fn denom(&self) -> u32 {
        let left_choice = (self.left_size - self.type_.size) as u32;
        let right_choice = (self.right_size - self.type_.size) as u32;
        binomial(left_choice, left_choice + right_choice)
    }
}

impl<F: Flag> Savable<Vec<CsMat<u32>>, F> for SplitCount<F> {
    fn filename(&self) -> String {
        format!(
            "split_{}_{}{}",
            self.left_size,
            self.right_size,
            self.type_.to_string_suffix()
        )
    }
    fn create(&self) -> Vec<CsMat<u32>> {
        let left = self.left_basis().get();
        let right = self.right_basis().get();
        let target = self.output_basis().get();
        count_split_tabulate(self.type_.size, &left, &right, &target)
    }
}

//================ Subflag count

/// .get() gives a matrix M where M\[i,j\] is the number of copies of i in j
#[derive(Clone, Debug)]
pub struct SubflagCount<F: Flag> {
    k: usize,
    n: usize,
    type_: Type<F>,
}

impl<F: Flag> SubflagCount<F> {
    pub fn make(k: usize, n: usize, type_: Type<F>) -> Self {
        assert!(type_.size <= k);
        assert!(k <= n);
        Self { k, n, type_ }
    }

    pub fn from_to(inner: Basis<F>, outer: Basis<F>) -> Self {
        assert_eq!(inner.t, outer.t);
        Self::make(inner.size, outer.size, inner.t)
    }

    pub fn inner_basis(&self) -> Basis<F> {
        Basis::make(self.k, self.type_)
    }

    fn outer_basis(&self) -> Basis<F> {
        Basis::make(self.n, self.type_)
    }

    pub fn denom(&self) -> u32 {
        let choices = (self.k - self.type_.size) as u32;
        let total = (self.n - self.type_.size) as u32;
        binomial(choices, total)
    }
}

impl<F: Flag> Savable<CsMat<u32>, F> for SubflagCount<F> {
    fn filename(&self) -> String {
        format!(
            "subflag_{}_to_{}{}",
            self.n,
            self.k,
            self.type_.to_string_suffix()
        )
    }
    fn create(&self) -> CsMat<u32> {
        let inner = self.inner_basis().get();
        let outer = self.outer_basis().get();
        count_subflag_tabulate(self.type_.size, &inner, &outer)
    }
}

// == unlabeling operators
/// Let F be the flag indexed by id on basis basis
/// this represents the unlabeling opearation that
/// sends the type `fully_typed(F)` to the flag `F`
#[derive(Debug, Clone)]
pub struct Unlabeling<F: Flag> {
    pub flag: usize,
    pub basis: Basis<F>,
}

impl<F: Flag> Unlabeling<F> {
    pub fn new(basis: Basis<F>, flag: usize) -> Self {
        Self { basis, flag }
    }
    pub fn total(t: Type<F>) -> Self {
        Self::new(Basis::new(t.size), t.id)
    }

    pub fn input_type(&self) -> Type<F>
    where
        F: Flag,
    {
        if self.basis.t == Type::empty() {
            Type::new(self.basis.size, self.flag) // !!! Do we assume something ?
        } else {
            let basis = self.basis.get();
            let unlab_basis = self.basis.without_type().get();
            let flag = &basis[self.flag];
            let unlab_id = unlab_basis.binary_search(&flag.canonical()).unwrap();
            Type::new(self.basis.size, unlab_id)
        }
    }

    pub fn output_type(&self) -> Type<F> {
        self.basis.t
    }
    /// Return the eta function of Razborov corresponding to
    /// the untyping operator.
    pub fn eta(&self) -> Vec<usize>
    where
        F: Flag,
    {
        let flag = &self.basis.get()[self.flag];
        let mut morphism = flag.morphism_to_canonical();
        morphism.truncate(self.basis.t.size);
        morphism
    }
}

#[derive(Debug, Clone)]
pub struct Unlabel<F: Flag> {
    pub unlabeling: Unlabeling<F>,
    pub size: usize,
}

impl<F: Flag> Savable<(Vec<usize>, Vec<u32>), F> for Unlabel<F> {
    fn filename(&self) -> String {
        format!(
            "unlabel_{}_id_{}_basis_{}{}",
            self.size,
            self.unlabeling.flag,
            self.unlabeling.basis.size,
            self.unlabeling.basis.t.to_string_suffix()
        )
    }
    fn create(&self) -> (Vec<usize>, Vec<u32>) {
        let in_basis = Basis::<F>::make(self.size, self.unlabeling.input_type()).get();
        let out_basis = Basis::<F>::make(self.size, self.unlabeling.output_type()).get();
        let eta = self.unlabeling.eta();
        (
            unlabeling_tabulate(&eta, &in_basis, &out_basis),
            unlabeling_count_tabulate(&eta, self.unlabeling.basis.size, &in_basis),
        )
    }
}

impl<F: Flag> Unlabel<F> {
    pub fn denom(&self) -> u32 {
        let new_type_size = self.unlabeling.basis.t.size;
        let old_type_size = self.unlabeling.basis.size;
        let choices = (old_type_size - new_type_size) as u32;
        let free_vertices = (self.size - new_type_size) as u32;
        product(free_vertices + 1 - choices, free_vertices)
    }
    pub fn output_basis(&self) -> Basis<F>
    where
        F: Flag,
    {
        Basis::new(self.size).with_type(self.unlabeling.output_type())
    }
    pub fn total(b: Basis<F>) -> Self {
        Self {
            unlabeling: Unlabeling::total(b.t),
            size: b.size,
        }
    }
}

#[derive(Debug, Clone)]
pub struct MulAndUnlabel<F: Flag> {
    pub split: SplitCount<F>,
    pub unlabeling: Unlabeling<F>,
}

impl<F: Flag> MulAndUnlabel<F> {
    pub fn invariant_classes(&self) -> InvariantClasses<F> {
        assert_eq!(self.split.left_size, self.split.right_size);
        let size = self.split.left_size;
        InvariantClasses(Unlabel {
            size,
            unlabeling: self.unlabeling,
        })
    }
    pub fn reduced(&self) -> ReducedByInvariant<F> {
        ReducedByInvariant(*self)
    }
}

impl<F: Flag> Display for MulAndUnlabel<F> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "Mul. and unlabel: {}x{}; {} -> {} (id {})",
            self.split.left_size,
            self.split.right_size,
            self.split.type_,
            self.unlabeling.basis.t,
            self.unlabeling.flag
        )
    }
}

impl<F: Flag> MulAndUnlabel<F> {
    fn unlabel(&self) -> Unlabel<F> {
        Unlabel {
            unlabeling: self.unlabeling,
            size: self.split.output_basis().size,
        }
    }
    pub fn output_basis(&self) -> Basis<F> {
        Basis::new(self.split.output_basis().size).with_type(self.unlabeling.output_type())
    }
    pub fn denom(&self) -> u32 {
        self.split.denom() * self.unlabel().denom()
    }
}

impl<F: Flag> Savable<Vec<CsMat<i64>>, F> for MulAndUnlabel<F> {
    fn filename(&self) -> String {
        format!(
            "{}_then_unlab_id_{}{}",
            self.split.filename(),
            self.unlabeling.flag,
            self.unlabeling.basis.t.to_string_suffix()
        )
    }
    fn create(&self) -> Vec<CsMat<i64>> {
        let (unlab_f, unlab_c) = self.unlabel().get();
        let mul = self.split.get();
        assert!(!mul.is_empty());
        let n = self.output_basis().get().len();
        let pre = pre_image(n, &unlab_f);
        let mut res = Vec::new();
        for pre_i in &pre {
            let mut res_i: CsMat<u32> = CsMat::zero(mul[0].shape());
            for &j in pre_i {
                let mut mul_j = mul[j].clone();
                mul_j *= unlab_c[j]; //
                res_i = &res_i + &mul_j; // can be optimized
            }
            res.push(res_i.map(|&v| v as i64))
        }
        debug_assert_eq!(res.len(), self.output_basis().get().len());
        debug_assert_eq!(
            res[0].shape(),
            (
                self.split.left_basis().get().len(),
                self.split.right_basis().get().len()
            )
        );
        res
    }
}

//
#[derive(Clone, Debug)]
pub struct InvariantClasses<F: Flag>(Unlabel<F>);

impl<F: Flag> Savable<Vec<usize>, F> for InvariantClasses<F> {
    fn filename(&self) -> String {
        format!(
            "invariant_classes_{}_id_{}_basis_{}{}",
            self.0.size,
            self.0.unlabeling.flag,
            self.0.unlabeling.basis.size,
            self.0.unlabeling.basis.t.to_string_suffix()
        )
    }
    fn create(&self) -> Vec<usize> {
        let unlabeling = self.0.unlabeling;
        let t = unlabeling.input_type();
        let flags = Basis::<F>::new(self.0.size).with_type(t).get();
        invariant_classes(&unlabeling.eta(), unlabeling.basis.size, &flags)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ReducedByInvariant<F: Flag>(MulAndUnlabel<F>);

impl<F: Flag> Savable<(Vec<CsMat<i64>>, Vec<CsMat<i64>>), F> for ReducedByInvariant<F> {
    fn filename(&self) -> String {
        format!("reduced_{}", self.0.filename())
    }
    fn create(&self) -> (Vec<CsMat<i64>>, Vec<CsMat<i64>>) {
        let class = self.0.invariant_classes().get();
        let (invariant_mat, antiinvariant_mat) = class_matrices(&class);
        let mul_and_unlabel = self.0.get();
        let mut res_inv = Vec::with_capacity(mul_and_unlabel.len());
        let mut res_anti = Vec::with_capacity(mul_and_unlabel.len());
        for m in mul_and_unlabel.into_iter() {
            let invariant = &(&invariant_mat.transpose_view() * &m) * &invariant_mat;
            let antiinvariant = if antiinvariant_mat.cols() == 0 {
                CsMat::zero((0, 0)) // avoiding a small bug of sprs
            } else {
                &(&antiinvariant_mat.transpose_view() * &m) * &antiinvariant_mat
            };
            res_inv.push(invariant);
            res_anti.push(antiinvariant);
        }
        (res_inv, res_anti)
    }
}

// Workaround to give Basis and Unlabeling the Copy trait
// (derive(Copy) does not derive the right bound when working
// with PhantomData)
impl<F: Flag> Copy for Type<F> {}
impl<F: Flag> Copy for Unlabeling<F> {}
impl<F: Flag> Copy for Basis<F> {}
impl<F: Flag> Copy for SplitCount<F> {}
impl<F: Flag> Copy for MulAndUnlabel<F> {}
// =======================

/// # Defining a quantum flags from a specified basis.
impl<F: Flag> Basis<F> {
    /// Sum of all flags of the basis.
    /// This is an expression of the 1 of the flag algebra.
    /// ```
    /// use flag_algebra::*;
    /// use flag_algebra::flags::Graph;
    ///
    /// let b = Basis::new(2);
    /// let one = b.one();
    /// let other: QFlag<i64, Graph> = b.random();
    /// assert_eq!(&one * &other, other.expand(Basis::new(4)));
    /// ```
    pub fn one<N>(self) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        assert!(F::HEREDITARY || self.size == self.t.size);
        let n = self.get().len();
        QFlag {
            basis: self,
            data: Array::from_elem(n, N::one()),
            scale: 1,
            expr: Expr::FromIndicator(|_, _| true, self),
        }
    }
    /// The zero vector in the specified basis.
    /// ```
    /// use flag_algebra::*;
    /// use flag_algebra::flags::Graph;
    ///
    /// let basis = Basis::new(3);
    /// let x: QFlag<i64, Graph> = basis.random();
    /// assert_eq!(basis.zero() + &x, x);
    /// ```
    pub fn zero<N>(self) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        QFlag {
            basis: self,
            data: Array::zeros(self.get().len()),
            scale: 1,
            expr: Expr::Zero,
        }
    }
    pub(crate) fn qflag_from_vec<N>(self, vec: Vec<N>) -> QFlag<N, F> {
        assert_eq!(self.get().len(), vec.len());
        QFlag {
            basis: self,
            data: Array::from(vec),
            scale: 1,
            expr: Expr::unknown(format!("from_vec({})", self.print_concise())),
        }
    }
    /// Return the formal sum of the flags of the basis `self`
    /// that satisfies some predicate `f`.
    ///
    /// The predicate `f` takes two arguments `g` and `sigma`, where `g` is a reference to
    /// the flag and `sigma` is the size of the labeled part.
    /// ```
    /// use flag_algebra::*;
    /// use flag_algebra::flags::Graph;
    ///
    /// // Sum of graphs of size 3 with an even number of edges
    /// let b = Basis::<Graph>::new(3);
    /// let sum = b.qflag_from_indicator(|g, _| g.edges().count() % 2 == 0 );
    ///
    /// let e3: QFlag<f64, Graph> = flag(&Graph::new(3, &[]));
    /// let p3 = flag(&Graph::new(3, &[(0, 1), (1, 2)]));
    /// assert_eq!(sum, e3 + &p3);
    ///
    /// /// Sum of the graphs of size 3 rooted on one vertex v
    /// /// where v has degree at least 1
    /// let t = Type::from_flag(&Graph::new(1, &[])); // Type for one vertex
    /// let basis = Basis::new(3).with_type(t);
    /// let sum: QFlag<f64, Graph> = basis.qflag_from_indicator(|g, _| g.edge(0, 1) || g.edge(0, 2) );
    /// ```
    pub fn qflag_from_indicator<N>(self, f: fn(&F, usize) -> bool) -> QFlag<N, F>
    where
        N: One + Zero,
    {
        let vec: Vec<_> = self
            .get()
            .iter()
            .map(|flag| {
                if f(flag, self.t.size) {
                    N::one()
                } else {
                    N::zero()
                }
            })
            .collect();
        QFlag {
            basis: self,
            data: Array::from(vec),
            scale: 1,
            expr: Expr::FromIndicator(f, self),
        }
    }
    /// Return the formal sum of `f(g)*g` on the flags `g` of the basis `self`.
    /// The second parameter of `f` is the size of the type of `g`.
    /// ```
    /// use flag_algebra::*;
    /// use flag_algebra::flags::Graph;
    ///
    /// // Sum of graphs of size 3 weighted by their number of edges
    /// let b = Basis::<Graph>::new(3);
    /// let sum: QFlag<f64, Graph>  = b.qflag_from_coeff(|g, _| g.edges().count() as f64 );
    /// ```

    pub fn qflag_from_coeff<N, M, P>(self, f: P) -> QFlag<N, F>
    where
        P: Fn(&F, usize) -> M + 'static,
        M: Into<N>,
    {
        let rc_f: Rc<dyn Fn(&F, usize) -> N> = Rc::new(move |a, b| f(a, b).into());
        self.qflag_from_coeff_rc(rc_f)
    }
    pub(crate) fn qflag_from_coeff_rc<N>(&self, f: Rc<dyn Fn(&F, usize) -> N>) -> QFlag<N, F> {
        let vec: Vec<_> = self.get().iter().map(|g| f(g, self.t.size)).collect();
        QFlag {
            basis: *self,
            data: Array::from(vec),
            scale: 1,
            expr: Expr::FromFunction(f, *self),
        }
    }
    pub fn random<N>(self) -> QFlag<N, F>
    where
        N: From<i16>,
    {
        let data: Vec<_> = (0..self.get().len())
            .map(|_| {
                let x: i16 = rand::random();
                N::from(x)
            })
            .collect();
        QFlag {
            basis: self,
            data: Array::from(data),
            scale: 1,
            expr: Expr::unknown(format!("random({})", self.print_concise())),
        }
    }
    pub(crate) fn flag_from_id<N>(self, id: usize) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        self.flag_from_id_with_base_size(id, self.get().len())
    }
    pub(crate) fn flag_from_id_with_base_size<N>(self, id: usize, size: usize) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        let mut res = QFlag {
            basis: self,
            data: Array::zeros(size),
            scale: 1,
            expr: Expr::Flag(id, self),
        };
        res.data[id] = N::one();
        res
    }
    /// Create a quantum flag containing exactly one flag.
    pub fn flag<N>(self, f: &F) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        assert_eq!(self.size, f.size());
        let flags = self.get();
        let mut data = Array::zeros(flags.len());
        let f1 = f.canonical_typed(self.t.size);
        let id = flags.binary_search(&f1).expect("Flag not found in basis");
        data[id] = N::one();
        QFlag {
            basis: self,
            data,
            scale: 1,
            expr: Expr::Flag(id, self),
        }
    }
    /// Returns the list of identifiers of all Square-and-unlabel operators
    /// that can be used in Cauchy-Schwarz inequalities for a problem on the basis `self`.
    pub fn all_cs(&self) -> Vec<MulAndUnlabel<F>> {
        let mut res = Vec::new();
        let n = self.size;
        // m: size of a cs basis
        for m in (n + self.t.size) / 2 + 1..=(2 * n - 1) / 2 {
            let sigma = 2 * m - n;
            let unlab_basis = Self::new(sigma).with_type(self.t);
            for unlab_id in 0..unlab_basis.get().len() {
                let unlabeling = Unlabeling::new(unlab_basis, unlab_id);
                let input_basis = Self::new(m).with_type(unlabeling.input_type());
                let split = SplitCount::from_input(&input_basis, &input_basis);
                res.push(MulAndUnlabel { split, unlabeling })
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flag::SubClass;
    use crate::flags::*;

    #[test]
    fn basis() {
        assert_eq!(Basis::<Graph>::new(5).get().len(), 34);
        assert_eq!(Basis::<Graph>::make(3, Type::new(1, 0)).get().len(), 6);
        assert_eq!(Basis::<Graph>::make(4, Type::new(2, 1)).get().len(), 20);
        //
        assert_eq!(Basis::<OrientedGraph>::new(3).get().len(), 7);
        assert_eq!(Basis::<OrientedGraph>::new(5).get().len(), 582);
        assert_eq!(
            Basis::<OrientedGraph>::make(3, Type::new(1, 0)).get().len(),
            15
        );
        assert_eq!(
            Basis::<OrientedGraph>::make(4, Type::new(2, 0)).get().len(),
            126
        );
        //
        assert_eq!(Basis::<DirectedGraph>::new(2).get().len(), 3);
        assert_eq!(Basis::<DirectedGraph>::new(3).get().len(), 16);
        //
        assert_eq!(
            Basis::<SubClass<OrientedGraph, TriangleFree>>::new(3)
                .get()
                .len(),
            6
        );
        assert_eq!(
            Basis::<SubClass<OrientedGraph, TriangleFree>>::new(5)
                .get()
                .len(),
            317
        );
        assert_eq!(
            Basis::<SubClass<OrientedGraph, TriangleFree>>::make(3, Type::new(2, 1))
                .get()
                .len(),
            8
        );
    }
    #[test]
    fn splitcount() {
        assert_eq!(56, SplitCount::<Graph>::make(5, 7, Type::new(2, 1)).denom());
        assert_eq!(3, SplitCount::<Graph>::make(3, 4, Type::new(2, 1)).denom());
        let _ = SplitCount::<Graph>::make(3, 2, Type::empty()).get();
        let _ = SplitCount::<Graph>::make(2, 3, Type::empty()).get();
        let _ = SplitCount::<Graph>::make(2, 3, Type::new(1, 0)).get();
        //
        let _ = SplitCount::<OrientedGraph>::make(2, 3, Type::new(1, 0)).get();
    }
    #[test]
    fn subflagcount() {
        assert_eq!(
            45,
            SplitCount::<Graph>::make(4, 10, Type::new(2, 0)).denom()
        );
        let _ = SubflagCount::<Graph>::make(2, 3, Type::new(1, 0)).get();
        let _ = SubflagCount::<Graph>::make(3, 4, Type::new(2, 1)).get();
        let _ = SubflagCount::<Graph>::make(5, 5, Type::empty()).get();
        let _ = SubflagCount::<Graph>::make(3, 5, Type::new(1, 0)).get();
    }
    #[test]
    fn unlabel() {
        let t = Type::new(3, 1);
        let unlabeling = Unlabeling::<Graph>::total(t);
        let size = 5;
        assert_eq!((Unlabel { unlabeling, size }).denom(), 60);
        let _ = (Unlabel { unlabeling, size }).get();
        //
        let b = Basis::new(3).with_type(Type::new(2, 1));
        let unlabeling = Unlabeling::<Graph>::new(b, 0);
        let _ = (Unlabel { unlabeling, size }).get();
    }
    #[test]
    fn mulandunlabeling() {
        let t = Type::new(2, 1);
        let unlabeling = Unlabeling::<Graph>::total(t);
        let split = SplitCount::make(3, 2, t);
        let _mau = (MulAndUnlabel { split, unlabeling }).get();
    }
    #[test]
    fn type_iterator() {
        assert_eq!(Type::<Graph>::types_with_size(4).count(), 11);
    }
    //     #[test]
    //     fn unlabeling_eta() {
    //         let b = Basis::<Graph>::new(5).with_type(Type::new(3, 1));
    //         let unlabeling = Unlabeling::new(b, 1);
    //         let eta = unlabeling.eta();
    //         let g = &unlabeling.basis.get()[unlabeling.flag];
    //         let t = unlabeling.basis.t;
    //         assert_eq!(g.induce(&eta), Basis::new(t.size).get()[t.id])
    //     }
}
