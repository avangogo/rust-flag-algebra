//! Computing and stocking operators of the flag algebra.

use crate::algebra::QFlag;
use crate::combinatorics::*;
use crate::density::*;
use crate::flag::Flag;
use crate::prettyprint::Expr;
use canonical_form::*;
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
use std::io::{stdout, BufReader, BufWriter, Write};
use std::marker::PhantomData;
use std::ops::*;
use std::path::*;

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
        filename.push(Path::new(&F::name()));
        filename.push(self.filename());
        let _ = filename.set_extension("dat");
        filename
    }
    /// (Re)create the object, save it in the corresponding file and return it.
    fn create_and_save(&self, path: &Path) -> A {
        println!("Creating {}", path.display());
        let value = self.create();
        let file = File::create(path).unwrap();
        let buf = BufWriter::new(file);
        bincode::serialize_into(buf, &value).unwrap();
        value
    }
    /// Load the object if the file exists and is valid.
    fn load(&self, path: &Path) -> Result<A, Box<Error>> {
        let file = File::open(&path)?;
        let mut buf = BufReader::new(file);
        let data = bincode::deserialize_from(&mut buf)?;
        Ok(data)
    }
    /// Function to automatically load the object if the file exists and
    /// is valid, or create and save it otherwise.
    fn get(&self) -> A {
        let path = self.file_path();
        if path.exists() {
            print!("Loading {}", path.display());
            stdout().flush().unwrap();
            match self.load(&path) {
                Ok(v) => {
                    println!(" done");
                    v
                }
                Err(e) => {
                    println!(" failed");
                    eprintln!("Error: {}", e);
                    self.create_and_save(&path)
                }
            }
        } else {
            let dir = path.parent().unwrap();
            match fs::create_dir_all(dir) {
                Ok(()) => {
                    self.create_and_save(&path)
                }
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
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct Type {
    /// Size of the type.
    pub size: usize,
    /// Index of the type in the list of unlabeled flags of this size.
    pub id: usize,
}

impl Type {
    /// Constructor for the type.
    pub fn new(size: usize, id: usize) -> Self {
        Type { size, id }
    }
    /// Create a type of size 0.
    pub fn empty() -> Self {
        Self::new(0, 0)
    }
    /// Write string that identify the type.
    pub fn to_string(&self) -> String {
        if *self == Self::empty() {
            String::new()
        } else {
            format!("type_{}_id_{}", self.size, self.id)
        }
    }
}

//============ Basis ===========

/// The set of flags of size size and type t
/// .get() returns an ordered vector containing all corresponding flags
#[derive(PartialEq, Eq, Debug, Hash)]
pub struct Basis<F> {
    /// Number of vertices of the flags of the basis.
    pub size: usize,
    /// Type of the flags of the basis.
    pub t: Type,
    phantom: PhantomData<F>,
}

impl<F: Flag> Basis<F> {
    /// Constructor for a basis.
    pub fn make(size: usize, t: Type) -> Self {
        assert!(t.size <= size);
        Basis {
            size,
            t,
            phantom: PhantomData,
        }
    }
    /// Basis of flag with `size` vertices and without type.
    pub fn new(size: usize) -> Self {
        Self::make(size, Type::empty())
    }
    /// Basis of flag with `size` vertices and same type as `self`.
    pub fn with_size(&self, size: usize) -> Self {
        Self::make(size, self.t)
    }
    /// Basis of flag with same size as `self` and type `t`.
    pub fn with_type(&self, t: Type) -> Self {
        Self::make(self.size, t)
    }
    /// Basis of flag with same size as `self` without type `t`.
    pub fn without_type(&self) -> Self {
        self.with_type(Type::empty())
    }
}

impl<F> Display for Basis<F> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        if self.t == Type::empty() {
            write!(f, "Flags of size {} without type", self.size)
        } else {
            write!(f, "Flags of size {} with type {:?})", self.size, self.t)
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
        format!("flags_{}_{}", self.size, self.t.to_string())
    }
    fn create(&self) -> Vec<F> {
        if self.t == Type::empty() {
            if self.size == 0 {
                F::all_flags(0)
            } else {
                F::generate_next(&self.with_size(self.size - 1).get())
            }
        } else {
            let type_basis = Basis::new(self.t.size).get();
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
#[derive(Clone, Copy, Debug)]
pub struct SplitCount<F> {
    left_size: usize,
    right_size: usize,
    type_: Type,
    phantom: PhantomData<F>,
}

impl<F: Flag> SplitCount<F> {
    pub fn make(left_size: usize, right_size: usize, type_: Type) -> Self {
        assert!(type_.size <= left_size);
        assert!(type_.size <= right_size);
        SplitCount {
            left_size,
            right_size,
            type_,
            phantom: PhantomData,
        }
    }

    pub fn from_input(left: &Basis<F>, right: &Basis<F>) -> Self {
        assert_eq!(left.t, right.t);
        Self::make(left.size, right.size, left.t)
    }

    pub fn left_basis(&self) -> Basis<F> {
        Basis::make(self.left_size, self.type_)
    }

    fn right_basis(&self) -> Basis<F> {
        Basis::make(self.right_size, self.type_)
    }

    fn output_basis(&self) -> Basis<F> {
        Basis::make(
            self.right_size + self.left_size - self.type_.size,
            self.type_,
        )
    }
    pub fn denom(&self) -> u64 {
        let left_choice = (self.left_size - self.type_.size) as u64;
        let right_choice = (self.right_size - self.type_.size) as u64;
        binomial(left_choice, left_choice + right_choice)
    }
}

impl<F: Flag> Savable<Vec<CsMat<u64>>, F> for SplitCount<F> {
    fn filename(&self) -> String {
        format!(
            "split_{}_{}_{}",
            self.left_size,
            self.right_size,
            self.type_.to_string()
        )
    }
    fn create(&self) -> Vec<CsMat<u64>> {
        let left = self.left_basis().get();
        let right = self.right_basis().get();
        let target = self.output_basis().get();
        count_split_tabulate(self.type_.size, &left, &right, &target)
    }
}

//================ Subflag count

/// .get() gives a matrix M where M[i,j] is the number of copies of i in j
#[derive(Clone, Copy, Debug)]
pub struct SubflagCount<F> {
    k: usize,
    n: usize,
    type_: Type,
    phantom: PhantomData<F>,
}

impl<F: Flag> SubflagCount<F> {
    pub fn make(k: usize, n: usize, type_: Type) -> Self {
        assert!(type_.size <= k && k <= n);
        SubflagCount {
            k,
            n,
            type_,
            phantom: PhantomData,
        }
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

    pub fn denom(&self) -> u64 {
        let choices = (self.k - self.type_.size) as u64;
        let total = (self.n - self.type_.size) as u64;
        binomial(choices, total)
    }
}

impl<F: Flag> Savable<CsMat<u64>, F> for SubflagCount<F> {
    fn filename(&self) -> String {
        format!(
            "subflag_{}_to_{}_{}",
            self.n,
            self.k,
            self.type_.to_string()
        )
    }
    fn create(&self) -> CsMat<u64> {
        let inner = self.inner_basis().get();
        let outer = self.outer_basis().get();
        count_subflag_tabulate(self.type_.size, &inner, &outer)
    }
}

// == unlabeling operators
/// Let F be the flag indexed by id on basis basis
/// this represents the unlabeling opearation that
/// sends the type fully_typed(F) to the flag F
#[derive(Debug)]
pub struct Unlabeling<F> {
    flag: usize,
    basis: Basis<F>,
}

impl<F: Flag> Unlabeling<F> {
    fn to_string(&self) -> String {
        format!(
            "unlab_id_{}_basis_{}_{}",
            self.flag,
            self.basis.size,
            self.basis.t.to_string()
        )
    }

    pub fn new(basis: Basis<F>, flag: usize) -> Self {
        Unlabeling { basis, flag }
    }

    pub fn total(t: Type) -> Self {
        Self::new(Basis::new(t.size), t.id)
    }

    pub fn input_type(&self) -> Type {
        if self.basis.t == Type::empty() {
            Type::new(self.basis.size, self.flag) // !!! Do we assume something ?
        } else {
            let basis = self.basis.get();
            let unlab_basis = self.basis.without_type().get();
            let flag = &basis[self.flag];
            let unlab_id = unlab_basis.binary_search(&canonical_form(flag)).unwrap();
            Type::new(self.basis.size, unlab_id)
        }
    }

    pub fn output_type(&self) -> Type {
        self.basis.t
    }

    pub fn eta(&self) -> Vec<usize> {
        let flag = &self.basis.get()[self.flag];
        let mut morphism = canonical_form_morphism(flag);
        morphism.resize(self.basis.t.size, 0);
        morphism
    }
}

#[derive(Debug, Clone, Copy)]
pub struct UnlabelingFlag<F> {
    pub unlabeling: Unlabeling<F>,
    pub size: usize,
}

impl<F: Flag> Savable<Vec<usize>, F> for UnlabelingFlag<F> {
    fn filename(&self) -> String {
        format!("flag_{}_{}", self.size, self.unlabeling.to_string())
    }
    fn create(&self) -> Vec<usize> {
        let in_basis = Basis::make(self.size, self.unlabeling.input_type()).get();
        let out_basis = Basis::make(self.size, self.unlabeling.output_type()).get();
        let eta = self.unlabeling.eta();
        unlabeling_tabulate::<F>(&eta, &in_basis, &out_basis)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MulAndUnlabeling<F> {
    split: SplitCount<F>,
    unlabeling: Unlabeling<F>,
}

impl<F> Display for MulAndUnlabeling<F> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "Mul. and unlabel: {}x{}; {:?} -> {:?} (id {})",
            self.split.left_size,
            self.split.right_size,
            self.split.type_,
            self.unlabeling.basis.t,
            self.unlabeling.flag
        )
    }
}

impl<F: Flag> MulAndUnlabeling<F> {
    pub fn new(split: SplitCount<F>, unlabeling: Unlabeling<F>) -> Self {
        debug_assert_eq!(split.type_, unlabeling.input_type());
        MulAndUnlabeling { split, unlabeling }
    }
    fn unlabeling_flag(&self) -> UnlabelingFlag<F> {
        UnlabelingFlag {
            unlabeling: self.unlabeling,
            size: self.split.output_basis().size,
        }
    }
    fn unlabeling_count(&self) -> UnlabelingCount<F> {
        UnlabelingCount {
            unlabeling: self.unlabeling,
            size: self.split.output_basis().size,
        }
    }
    pub fn output_basis(&self) -> Basis<F> {
        Basis::new(self.split.output_basis().size).with_type(self.unlabeling.output_type())
    }
    pub fn denom(&self) -> u64 {
        self.split.denom() * self.unlabeling_count().denom()
    }
}

impl<F: Flag> Savable<Vec<CsMat<u64>>, F> for MulAndUnlabeling<F> {
    fn filename(&self) -> String {
        format!(
            "{}_then_{}",
            self.split.filename(),
            self.unlabeling.to_string()
        )
    }
    fn create(&self) -> Vec<CsMat<u64>> {
        // can be optimized
        let unlab_c = self.unlabeling_count().get();
        let unlab_f = self.unlabeling_flag().get();
        let mul = self.split.get();
        assert!(!mul.is_empty());
        let n = self.output_basis().get().len();
        let pre = pre_image(n, &unlab_f);
        let mut res = Vec::new();
        for pre_i in pre.iter() {
            let mut res_i = CsMat::zero(mul[0].shape());
            for &j in pre_i {
                res_i = &res_i + &(&mul[j] * unlab_c[j]); // can be optimized
            }
            res.push(res_i)
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

// Workaround to give Basis the Copy trait
// (derive(Copy) does not to work well with PhantomData)
impl<F> Clone for Basis<F> {
    fn clone(&self) -> Self {
        Basis {
            size: self.size,
            t: self.t,
            phantom: PhantomData,
        }
    }
}
impl<F> Clone for Unlabeling<F> {
    fn clone(&self) -> Self {
        Unlabeling {
            flag: self.flag,
            basis: self.basis,
        }
    }
}
impl<F> Copy for Basis<F> {}
impl<F> Copy for Unlabeling<F> {}

// ==============
#[derive(Clone, Copy, Debug)]
pub struct UnlabelingCount<F> {
    pub unlabeling: Unlabeling<F>,
    pub size: usize,
}

impl<F: Flag> Savable<Vec<u64>, F> for UnlabelingCount<F> {
    fn filename(&self) -> String {
        format!("count_{}_{}", self.size, self.unlabeling.to_string())
    }
    fn create(&self) -> Vec<u64> {
        let in_basis: Vec<F> = Basis::make(self.size, self.unlabeling.input_type()).get();
        unlabeling_count_tabulate(
            &self.unlabeling.eta(),
            self.unlabeling.basis.size,
            &in_basis,
        )
    }
}

impl<F> UnlabelingCount<F> {
    pub fn denom(&self) -> u64 {
        let initial_type_size = self.unlabeling.basis.t.size;
        let choices = (self.unlabeling.basis.size - initial_type_size) as u64;
        let total = (self.size - initial_type_size) as u64;
        product(total, total - choices - 1)
    }
}

// =======================
///Constructing quantum graphs on a basis
// FIXME : To be removed ?
impl<F> Basis<F>
where
    F: Flag,
{
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

    pub fn one<N>(self) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        let n = self.get().len();
        QFlag {
            basis: self,
            data: Array::from_elem(n, N::one()),
            scale: 1,
            expr: Expr::One,
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
            data: Array::from_vec(data),
            scale: 1,
            expr: Expr::Num(String::from("random")),
        }
    }

    pub fn flag_from_id<N>(self, id: usize) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        self.flag_from_id_with_base_size(id, self.get().len())
    }

    pub fn flag_from_id_with_base_size<N>(self, id: usize, size: usize) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        let mut res = QFlag {
            basis: self,
            data: Array::zeros(size),
            scale: 1,
            expr: Expr::Num(String::from("from_id")),
        };
        res.data[id] = N::one();
        res
    }

    pub fn flag<N>(self, f: &F) -> QFlag<N, F>
    where
        N: Num + Clone,
    {
        assert_eq!(self.size, f.size());
        let flags = self.get();
        let mut data = Array::zeros(flags.len());
        let f1 = canonical_form_typed(f, self.t.size);
        data[flags.binary_search(&f1).expect("Flag not found in basis")] = N::one();
        QFlag {
            basis: self,
            data,
            scale: 1,
            expr: Expr::Num(String::from("flag")),
        }
    }
    pub fn from_vec<N>(self, vec: Vec<N>) -> QFlag<N, F> {
        assert_eq!(self.get().len(), vec.len());
        QFlag {
            basis: self,
            data: Array::from_vec(vec),
            scale: 1,
            expr: Expr::Num(String::from("from_vec")),
        }
    }
    pub fn from_indicator<M, N, P>(self, mut f: P) -> QFlag<N, F>
    where
        M: Into<N>,
        P: FnMut(&F, usize) -> M,
    {
        let flags = self.get();
        let mut vec = Vec::new();
        for g in flags.iter() {
            vec.push(f(g, self.t.size).into())
        }
        QFlag {
            basis: self,
            data: Array::from_vec(vec),
            scale: 1,
            expr: Expr::Num(String::from("sum f(F)F")),
        }
    }
    pub fn all_cs(&self) -> Vec<MulAndUnlabeling<F>> {
        let mut res = Vec::new();
        let n = self.size;
        // m: size of a cs basis
        for m in (n + self.t.size + 1) / 2..=(2 * n - 1) / 2 {
            let sigma = 2 * m - n;
            let unlab_basis = Basis::<F>::new(sigma).with_type(self.t);
            for unlab_id in 0..unlab_basis.get().len() {
                let unlab = Unlabeling::new(unlab_basis, unlab_id);
                let input_basis = Basis::new(m).with_type(unlab.input_type());
                let split = SplitCount::from_input(&input_basis, &input_basis);
                res.push(MulAndUnlabeling::new(split, unlab))
            }
        }
        res
    }
}

pub fn flag_typed<N, F>(f: &F, type_size: usize) -> QFlag<N, F>
where
    N: Num + Clone,
    F: Flag,
{
    let flag = canonical_form_typed(f, type_size);
    let type_ = flag.induce(&(0..type_size).collect::<Vec<_>>()); // type
    let type_basis = Basis::new(type_size);
    let type_id = type_basis
        .get()
        .binary_search(&type_)
        .expect("Flag not found in basis");
    let t = Type::new(type_size, type_id);
    let basis = Basis::new(f.size()).with_type(t);
    basis.flag(&flag)
}

pub fn flag<N, F>(f: &F) -> QFlag<N, F>
where
    N: Num + Clone,
    F: Flag,
{
    flag_typed(f, 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::*;

    #[test]
    fn basis() {
        assert_eq!(Basis::<Graph>::new(5).get().len(), 34);
        assert_eq!(Basis::<Graph>::make(3, Type::new(1, 0)).get().len(), 6);
        assert_eq!(Basis::<Graph>::make(4, Type::new(2, 1)).get().len(), 20);
        //
        assert_eq!(Basis::<Digraph>::new(3).get().len(), 7);
        assert_eq!(Basis::<Digraph>::new(5).get().len(), 582);
        assert_eq!(Basis::<Digraph>::make(3, Type::new(1, 0)).get().len(), 15);
        assert_eq!(Basis::<Digraph>::make(4, Type::new(2, 0)).get().len(), 126);
    }
    #[test]
    fn splitcount() {
        assert_eq!(56, SplitCount::<Graph>::make(5, 7, Type::new(2, 1)).denom());
        let _ = SplitCount::<Graph>::make(3, 2, Type::empty()).get();
        let _ = SplitCount::<Graph>::make(2, 3, Type::empty()).get();
        let _ = SplitCount::<Graph>::make(2, 3, Type::new(1, 0)).get();
        //
        let _ = SplitCount::<Digraph>::make(2, 3, Type::new(1, 0)).get();
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
    fn unlabeling() {
        let t = Type { size: 3, id: 1 };
        let unlabeling = Unlabeling::<Graph>::total(t);
        let size = 5;
        let _ = (UnlabelingCount { unlabeling, size }).get();
        let _ = (UnlabelingFlag { unlabeling, size }).get();
        //
        let b = Basis::new(3).with_type(Type { size: 2, id: 1 });
        let unlabeling = Unlabeling::<Graph>::new(b, 0);
        let _ = (UnlabelingCount { unlabeling, size }).get();
        let _ = (UnlabelingFlag { unlabeling, size }).get();
    }
    #[test]
    fn mulandunlabeling() {
        let t = Type { size: 2, id: 1 };
        let unlabeling = Unlabeling::<Graph>::total(t);
        let _mau = MulAndUnlabeling::new(SplitCount::make(3, 2, t), unlabeling).get();
    }
    // #[test]
    // fn unlabeling_eta() {
    //     let b = Basis::<Graph>::new(3).with_type(Type::new(2,1));
    //     let unlabeling = Unlabeling::new(b, 1);
    //     let u = unlabeling.eta();
    // }
}
