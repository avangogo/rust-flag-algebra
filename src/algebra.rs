//! Manipulation of vectors and inequalities in a flag algebra.

use crate::expr::{Expr, Names, VarRange};
use crate::flag::Flag;
use crate::operator::*;
use ndarray::{Array1, ScalarOperand};
use num::{pow::Pow, FromPrimitive, Integer, Num};
use sprs::{CsMat, CsMatView, CsVec, MulAcc, TriMat};
use std::fmt::*;
use std::ops::*;

/// An element of a flag algebra.
#[derive(Clone, Debug)]
pub struct QFlag<N, F: Flag> {
    /// Basis of the space where the element lives. This corresponds to the size and type of the flags.
    pub basis: Basis<F>,
    /// The vector of the element in the corresponding basis is `(1/self.scale).self.data`.
    pub data: Array1<N>,
    /// Scaling factor of the vector.
    pub scale: u64,
    /// Expression recording how the vector was computed.
    pub expr: Expr<N, F>,
}

// equality for QFlags
impl<N, F> PartialEq for QFlag<N, F>
where
    N: Num + FromPrimitive + Clone,
    F: Flag,
{
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.basis, other.basis);
        assert_eq!(self.data.len(), other.data.len());
        //
        let s1 = N::from_u64(self.scale).unwrap();
        let s2 = N::from_u64(other.scale).unwrap();
        self.data
            .iter()
            .zip(other.data.iter())
            .all(|(x, y)| x.clone() * s2.clone() == y.clone() * s1.clone())
    }
}

// ==================== operarions on Qflags ===========

/// Arithmetic to to put two scaled vectors on same denominator
///
/// If `f1 = v1 / scale1` and `f2 = v2 / scale2`, then
/// `matching_scales(scale1, scale2)` returns `(c1, c2, scale)`
/// such that  `f1 = v1 * c1 / scale` and `f2 = v2 * c2 / scale`.
fn matching_scales<N>(scale1: u64, scale2: u64) -> (N, N, u64)
where
    N: FromPrimitive,
{
    let gcd = scale1.gcd(&scale2);
    let c1 = N::from_u64(scale2 / gcd).unwrap();
    let c2 = N::from_u64(scale1 / gcd).unwrap();
    let scale = (scale1 / gcd) * scale2;
    (c1, c2, scale)
}

impl<N, F> Add<&Self> for QFlag<N, F>
where
    N: Clone + FromPrimitive + Num + ScalarOperand,
    F: Flag,
{
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        assert_eq!(self.basis, other.basis);
        assert_eq!(self.data.len(), other.data.len());
        let (a1, a2, scale) = matching_scales::<N>(self.scale, other.scale);
        QFlag {
            basis: self.basis,
            data: self.data * a1 + &other.data * a2,
            scale,
            expr: self.expr + other.expr.clone(),
        }
    }
}

impl<N, F> Add for QFlag<N, F>
where
    N: Clone + Num + FromPrimitive + ScalarOperand,
    F: Flag,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        self + &other
    }
}

impl<'a, N, F> Sub for &'a QFlag<N, F>
where
    N: Clone + Num + FromPrimitive + ScalarOperand,
    F: Flag,
{
    type Output = QFlag<N, F>;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(self.basis, other.basis);
        assert_eq!(self.data.len(), other.data.len());
        let (a1, a2, scale) = matching_scales::<N>(self.scale, other.scale);
        QFlag {
            basis: self.basis,
            data: &self.data * a1 - &other.data * a2,
            scale,
            expr: self.expr.clone() - other.expr.clone(),
        }
    }
}

impl<N, F> Sub for QFlag<N, F>
where
    N: Clone + Num + FromPrimitive + ScalarOperand,
    F: Flag,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        &self - &other
    }
}

impl<N, F: Flag> Neg for QFlag<N, F>
where
    N: Clone + Neg<Output = N>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            basis: self.basis,
            data: -self.data,
            scale: self.scale,
            expr: -self.expr,
        }
    }
}

impl<'a, N, F> Neg for &'a QFlag<N, F>
where
    N: Clone + Neg<Output = N>,
    F: Flag,
{
    type Output = QFlag<N, F>;

    fn neg(self) -> Self::Output {
        QFlag {
            basis: self.basis,
            data: -self.data.clone(),
            scale: self.scale,
            expr: -self.expr.clone(),
        }
    }
}

// Right scalar multiplication (it is not possible to implement it on left)
impl<N, F> Mul<N> for QFlag<N, F>
where
    N: Num + ScalarOperand + Display,
    F: Flag,
{
    type Output = Self;

    fn mul(self, rhs: N) -> Self::Output {
        Self {
            expr: Expr::num(&rhs) * self.expr.clone(),
            basis: self.basis,
            data: self.data * rhs,
            scale: self.scale,
        }
    }
}

impl<N, F> Pow<usize> for &QFlag<N, F>
where
    N: Num + Clone + FromPrimitive + Display,
    F: Flag,
{
    type Output = QFlag<N, F>;

    fn pow(self, n: usize) -> QFlag<N, F> {
        match n {
            0 => self.basis.with_size(self.basis.t.size).one(),
            1 => self.clone(),
            n => {
                let mut res = self * self;
                for _ in 2..n {
                    res = &res * self
                }
                res
            }
        }
    }
}

impl<N, F: Flag> Display for IneqMeta<N, F>
where
    N: Display,
{
    fn fmt(&self, f: &mut Formatter) -> Result {
        write!(
            f,
            "{}\t{} {}",
            self.flag_expr,
            if self.equality { '=' } else { '≥' },
            self.bound_expr
        )
    }
}

impl<N, F: Flag> Display for Ineq<N, F>
where
    N: Display,
{
    fn fmt(&self, f: &mut Formatter) -> Result {
        self.meta.fmt(f)
    }
}

// =================
fn quadratic_form<N>(lhs: &Array1<N>, matrix: &CsMat<u32>, rhs: &Array1<N>) -> N
where
    N: Num + Clone + FromPrimitive,
{
    assert_eq!(lhs.len(), matrix.rows());
    assert_eq!(rhs.len(), matrix.cols());
    let mut res = N::zero();
    for (v, (i, j)) in matrix {
        res = res + (N::from_u32(*v).unwrap() * lhs[i].clone() * rhs[j].clone());
    }
    res
}

fn vector_matrix_mul<N>(matrix: &CsMatView<u32>, vec: &Array1<N>) -> Array1<N>
where
    N: Num + Clone + FromPrimitive,
{
    assert_eq!(vec.len(), matrix.cols());
    let mut res: Array1<N> = Array1::zeros(matrix.rows());
    for (&v, (i, j)) in matrix {
        res[i] = res[i].clone() + N::from_u32(v).unwrap() * vec[j].clone();
    }
    res
}

fn multiply<N>(lhs: &Array1<N>, table: &[CsMat<u32>], rhs: &Array1<N>) -> Array1<N>
where
    N: Num + Clone + FromPrimitive,
{
    let mut res = Array1::<N>::zeros(table.len());
    for (i, matrix) in table.iter().enumerate() {
        res[i] = quadratic_form(lhs, matrix, rhs);
    }
    res
}

// Conversions between dense and sparse arrays
// Dense to sparse
fn csvec_from_array<N>(array: &Array1<N>) -> CsVec<N>
where
    N: Num + Clone,
{
    let mut res = CsVec::empty(array.len());
    for (i, val) in array.iter().enumerate() {
        if val != &N::zero() {
            res.append(i, val.clone())
        }
    }
    res
}

// Sparse to dense
fn array_from_csvec<N>(csvec: &CsVec<N>) -> Array1<N>
where
    N: Num + Clone,
{
    let mut res = vec![N::zero(); csvec.dim()];
    csvec.scatter(&mut res);
    Array1::from(res)
}

/// Flag operator function where the data from the flag algebra is given in input
impl<N, F> QFlag<N, F>
where
    N: Num + Clone + FromPrimitive,
    F: Flag,
{
    fn raw_expand(&self, operator: &CsMat<u32>, outbasis: Basis<F>, denom: u32) -> Self {
        Self {
            basis: outbasis,
            data: vector_matrix_mul(&operator.view(), &self.data),
            scale: self.scale * denom as u64,
            expr: self.expr.clone(),
        }
    }
    fn raw_multiply(&self, table: &[CsMat<u32>], other: &Self, denom: u32) -> Self {
        assert_eq!(self.basis.t, other.basis.t);
        Self {
            basis: self.basis * other.basis,
            data: multiply(&self.data, table, &other.data),
            scale: self.scale * denom as u64 * other.scale,
            expr: Expr::mul(self.expr.clone(), other.expr.clone()),
        }
    }
    fn raw_untype(
        &self,
        untype_flag: &[usize],
        untype_count: &[u32],
        outbasis: Basis<F>,
        outbasis_size: usize,
        denom: u32,
    ) -> Self {
        assert_eq!(untype_flag.len(), untype_count.len());
        let mut data = Array1::<N>::zeros(outbasis_size);
        for (i, v) in self.data.iter().enumerate() {
            data[untype_flag[i]] =
                data[untype_flag[i]].clone() + v.clone() * N::from_u32(untype_count[i]).unwrap()
        }
        Self {
            basis: outbasis,
            data,
            scale: self.scale * denom as u64,
            expr: self.expr.clone().unlab(),
        }
    }
}

fn untype_matrix<N>(untype_flag: &[usize], untype_count: &[u32], outbasis_size: usize) -> CsMat<N>
where
    N: Num + FromPrimitive + Clone,
{
    let inbasis_size = untype_flag.len();
    let shape = (outbasis_size, inbasis_size);
    let mut trimat = TriMat::with_capacity(shape, inbasis_size);
    for i in 0..untype_flag.len() {
        trimat.add_triplet(untype_flag[i], i, N::from_u32(untype_count[i]).unwrap())
    }
    trimat.to_csr()
}

/// # Flag algebra operations

impl<'a, N, F> Mul for &'a QFlag<N, F>
where
    N: Num + Clone + FromPrimitive,
    F: Flag,
{
    type Output = QFlag<N, F>;

    fn mul(self, other: Self) -> QFlag<N, F> {
        let split = SplitCount::from_input(&self.basis, &other.basis);
        self.raw_multiply(&split.get(), other, split.denom())
    }
}

impl<N, F> Mul for QFlag<N, F>
where
    N: Num + Clone + FromPrimitive,
    F: Flag,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        &self * &other
    }
}

impl<N, F: Flag> QFlag<N, F> {
    /// Projection to a basis of larger flag.
    pub fn expand(&self, outbasis: Basis<F>) -> Self
    where
        N: Num + Clone + FromPrimitive,
    {
        let subflag = SubflagCount::from_to(self.basis, outbasis);
        self.raw_expand(&subflag.get(), outbasis, subflag.denom())
    }
    /// Unlabeling operator 〚.〛 to the flag algebra of completly unlabeled flags.
    pub fn untype(&self) -> Self
    where
        N: Num + Clone + FromPrimitive,
    {
        let unlabeling = Unlabeling::<F>::total(self.basis.t);
        let size = self.basis.size;
        let outbasis = self.basis.with_type(Type::empty());
        let unlabel = Unlabel { unlabeling, size };
        let (unlab_flag, unlab_count) = unlabel.get();
        self.raw_untype(
            &unlab_flag,
            &unlab_count,
            outbasis,
            outbasis.get().len(),
            unlabel.denom(),
        )
    }
    /// Return the same element with modified pretty-print expression
    pub fn with_expr(mut self, expr: Expr<N, F>) -> Self {
        self.expr = expr;
        self
    }
    /// Name the vector for the purpose of pretty-printing
    pub fn named(mut self, name: String) -> Self {
        self.expr = self.expr.named(name);
        self
    }
    /// Divide the elements of the vector by the scale and set the scale to 1
    pub fn no_scale(mut self) -> Self
    where
        N: FromPrimitive + DivAssign<N> + ScalarOperand,
    {
        self.data /= N::from_u64(self.scale).unwrap();
        self.scale = 1;
        self
    }
    pub fn map<G, M>(&self, g: G) -> QFlag<M, F>
    where
        G: Fn(&N) -> M,
    {
        QFlag {
            basis: self.basis,
            data: self.data.map(&g),
            scale: self.scale,
            expr: self.expr.map(&g),
        }
    }
}

/// # Creating inequalies from `QFlag`
impl<N, F> QFlag<N, F>
where
    N: Num + FromPrimitive + Clone + Display,
    F: Flag,
{
    /// Return the inequality "`self` ≥ `x`".
    pub fn at_least(&self, x: N) -> Ineq<N, F> {
        Ineq {
            meta: IneqMeta {
                basis: self.basis,
                flag_expr: self.expr.clone(),
                bound_expr: Expr::num(&x),
                equality: false,
                forall: None,
                scale: self.scale,
            },
            data: vec![IneqData {
                flag: csvec_from_array(&self.data),
                bound: x * N::from_u64(self.scale).unwrap(),
            }],
        }
    }
    /// Return the inequality "`self` ≤ `x`".
    pub fn at_most(&self, x: N) -> Ineq<N, F>
    where
        N: Clone + Neg<Output = N>,
    {
        (-self.clone()).at_least(-x)
    }
    /// Return the inequality "`self` ≥ `0`".
    pub fn non_negative(&self) -> Ineq<N, F>
    where
        N: Num,
    {
        self.at_least(N::zero())
    }
    /// Return the equality "`self` = `n`".
    pub fn equal(self, n: N) -> Ineq<N, F>
    where
        N: Clone + Neg<Output = N>,
    {
        self.at_least(n).equality()
    }
}

/// Return the inequalities expressing that the sum of the flags of `basis`
/// is equal to one.
pub fn total_sum_is_one<N, F>(basis: Basis<F>) -> Ineq<N, F>
where
    F: Flag,
    N: Num + Clone + Neg<Output = N> + FromPrimitive + Display,
{
    basis.one().equal(N::one())
}

/// Return the inequalities expressing that the flags of `basis`
/// are larger than zero.
pub fn flags_are_nonnegative<N, F>(basis: Basis<F>) -> Ineq<N, F>
where
    F: Flag,
    N: Num + Clone + Neg<Output = N>,
{
    let n = basis.get().len();
    let mut data = Vec::with_capacity(n);
    for i in 0..n {
        let mut flag = CsVec::empty(n);
        flag.append(i, N::one());
        data.push(IneqData {
            flag,
            bound: N::zero(),
        })
    }
    let meta = IneqMeta {
        basis,
        flag_expr: Expr::Var(0).named(format!("flag(:{})", basis.print_concise())),
        bound_expr: Expr::Zero,
        equality: false,
        forall: Some(VarRange::InBasis(basis)),
        scale: 1,
    };
    Ineq { meta, data }
}

//============
#[derive(Clone, Debug)]
/// Contains informations about a set of inequalities of a flag algebra.
pub(crate) struct IneqMeta<N, F: Flag> {
    /// Basis in which the inequality is expressed.
    /// This correspond to the type and size of the flags.
    pub basis: Basis<F>,
    /// Expression recording how the left sides of the inequalities where constructed.
    pub flag_expr: Expr<N, F>,
    ///
    forall: Option<VarRange<F>>,
    /// Expression recording how the right sides where constructed.
    pub bound_expr: Expr<N, F>,
    /// True if the inequality is an equality
    pub equality: bool,
    ///
    scale: u64,
}

impl<N: Clone, F: Flag> IneqMeta<N, F> {
    fn opposite(self) -> Self {
        Self {
            basis: self.basis,
            flag_expr: self.flag_expr.neg(),
            bound_expr: self.bound_expr.neg(),
            equality: self.equality,
            forall: self.forall,
            scale: self.scale,
        }
    }
    fn one_sided_expr(&self) -> Expr<N, F> {
        Expr::sub(self.flag_expr.clone(), self.bound_expr.clone())
    }

    fn multiply(&self, rhs_basis: Basis<F>, rhs_expr: Expr<N, F>) -> Self {
        let forall = if let Expr::Var(_) = rhs_expr {
            match self.forall {
                None => Some(VarRange::InBasis(rhs_basis)),
                Some(_) => unimplemented!(),
            }
        } else {
            self.forall.clone()
        };
        Self {
            basis: self.basis * rhs_basis,
            flag_expr: Expr::mul(self.one_sided_expr(), rhs_expr),
            bound_expr: Expr::Zero,
            equality: self.equality,
            forall,
            scale: self.scale * SplitCount::from_input(&self.basis, &rhs_basis).denom() as u64,
        }
    }

    fn untype(&self) -> Self {
        Self {
            basis: self.basis.with_type(Type::empty()),
            flag_expr: Expr::unlab(self.flag_expr.clone()),
            bound_expr: self.bound_expr.clone(),
            equality: self.equality,
            forall: self.forall.clone(),
            scale: self.scale * Unlabel::total(self.basis).denom() as u64,
        }
    }

    pub(crate) fn latex(&self, names: &mut Names<N, F>) -> String
    where
        N: Display,
    {
        format!(
            "{}{} {} {}",
            if let Some(ref range) = self.forall {
                range.latex(names)
            } else {
                String::new()
            },
            self.flag_expr.latex(names),
            if self.equality { "=" } else { "\\geq" },
            self.bound_expr.latex(names),
        )
    }
}

#[derive(Clone, Debug)]
/// Contains the vector and the bound of one inequality in a flag algebra.
/// This inequality has the form `self.flag  ≥ self.bound`.
/// Expression recording how the left sides where constructed.
pub(crate) struct IneqData<N> {
    /// Vector of the left side in the corresponding flag basis.
    pub flag: CsVec<N>,
    /// Number on the right side of the inequality.
    pub bound: N,
}

impl<N> IneqData<N>
where
    N: Num + Clone,
{
    fn opposite(self) -> Self
    where
        N: Neg<Output = N>,
    {
        let mut flag = self.flag;
        flag.map_inplace(|x| -x.clone());
        Self {
            flag,
            bound: -self.bound,
        }
    }
    fn one_sided(self) -> Self
    where
        N: Neg<Output = N>,
    {
        if self.bound == N::zero() {
            self
        } else {
            // Can be slow
            // compute self.flag - self.bound
            let n = self.flag.dim();
            let mut flag = CsVec::empty(n);
            flag.reserve(n);
            let mut next_j = 0;
            for (i, val) in self.flag.iter() {
                for j in next_j..i {
                    flag.append(j, -self.bound.clone())
                }
                flag.append(i, val.clone() - self.bound.clone());
                next_j = i + 1;
            }
            for j in next_j..n {
                flag.append(j, -self.bound.clone())
            }
            Self {
                flag,
                bound: N::zero(),
            }
        }
    }
    fn untype(&self, untype_matrix: &CsMat<N>, denom: u32) -> Self
    where
        N: Clone + Num + Default + MulAcc + Send + Sync + FromPrimitive,
    {
        Self {
            flag: untype_matrix * &self.flag,
            bound: self.bound.clone() * N::from_u32(denom).unwrap(),
        }
    }
    // From profiling: Memory allocation here could be optimized
    fn multiply_by_all(self, table: &[CsMat<N>], acc: &mut Vec<Self>)
    where
        N: Num + Clone + Send + Sync + MulAcc + Default + Neg<Output = N>,
    {
        if let Some(other_size) = table.first().map(|mat| mat.cols()) {
            let one_sided = self.one_sided();
            let mut flags: Vec<CsVec<N>> = vec![CsVec::empty(table.len()); other_size];
            for (i, mat) in table.iter().enumerate() {
                let vec: CsVec<N> = &mat.transpose_view() * &one_sided.flag.view();
                for (j, val) in vec.iter() {
                    flags[j].append(i, val.clone())
                }
            }
            for flag in flags {
                let ineq_data = Self {
                    flag,
                    bound: N::zero(),
                };
                acc.push(ineq_data)
            }
        }
    }
}

#[derive(Clone, Debug)]
/// A set of bounds on elements of a flag algebra.
///
/// This correpond to a set of inequalities constructed in a similar way.
pub struct Ineq<N, F: Flag> {
    /// Common information about the set of inequalities.
    pub(crate) meta: IneqMeta<N, F>,
    /// List of data of the inequalities in the set.
    pub(crate) data: Vec<IneqData<N>>,
}

impl<N, F> Ineq<N, F>
where
    N: Clone + Num,
    F: Flag,
{
    /// If self is "`f ≥ x`", returns "`f ≤ x`".
    pub fn opposite(self) -> Self
    where
        N: Neg<Output = N>,
    {
        Self {
            meta: self.meta.opposite(),
            data: self.data.into_iter().map(|x| x.opposite()).collect(),
        }
    }
    /// If self is "`f ≥ x`", returns "`f = x`".
    pub fn equality(mut self) -> Self {
        self.meta.equality = true;
        self
    }

    /// If self is "`f ≥ x`", returns "`f ≥ x - eps`".
    pub fn relaxed(mut self, eps: N) -> Self
    where
        N: SubAssign,
    {
        for ineq in &mut self.data {
            ineq.bound -= eps.clone()
        }
        self
    }
    /// Return the flag member of the ith inequality
    pub fn lhs(&self, i: usize) -> QFlag<N, F> {
        assert!(i < self.data.len());
        QFlag {
            basis: self.meta.basis,
            data: array_from_csvec(&self.data[i].flag),
            scale: self.meta.scale,
            expr: self.meta.flag_expr.substitute_option(&self.meta.forall, i),
        }
    }
    pub fn check(&self)
    where
        N: Debug + Neg<Output = N> + Clone + FromPrimitive + ScalarOperand + Display,
    {
        for i in 0..self.data.len() {
            let x = self.lhs(i);
            assert_eq!(x, x.expr.eval());
        }
    }
}

impl<N, F> Ineq<N, F>
where
    N: Num + Clone + Send + Sync + Default + FromPrimitive + AddAssign + std::iter::Sum + MulAcc,
    F: Flag,
{
    /// If self is "`f` ≥ `x`", return the projection "`〚f〛 ≥ x`".
    pub fn untype(&self) -> Self {
        let unlabeling = Unlabeling::<F>::total(self.meta.basis.t);
        let size = self.meta.basis.size;
        let unlabel = Unlabel { unlabeling, size };
        let (unlab_f, unlab_c) = unlabel.get();
        let outbasis_size = unlabel.output_basis().get().len();
        let unlab_matrix = untype_matrix(&unlab_f, &unlab_c, outbasis_size);
        let denom = unlabel.denom();
        //
        let mut data = Vec::new();
        for i in &self.data {
            let f = i.untype(&unlab_matrix, denom);
            data.push(f)
        }
        Self {
            meta: self.meta.untype(),
            data,
        }
    }
    /// If self is "`f` ≥ `x`", return the set of inequalities "`f*g ≥ x.g`",
    /// where `g` is chosen such that `f*g` is a vector of `outbasis`.
    pub fn multiply_by_all(self, outbasis: Basis<F>) -> Self
    where
        N: Neg<Output = N>,
    {
        let b = outbasis / self.meta.basis;
        let splitcount = SplitCount::from_input(&self.meta.basis, &b);
        let table: Vec<CsMat<N>> = splitcount
            .get()
            .iter()
            .map(|m| m.map(|&x| N::from_u32(x).unwrap()))
            .collect();
        //
        let mut data = Vec::new();
        for ineq in self.data {
            ineq.multiply_by_all(&table, &mut data)
        }
        //
        Self {
            data,
            meta: self.meta.multiply(b, Expr::Var(0)),
        }
    }
    /// If self is "`f` ≥ `x`", return the set of inequalities "`〚f*g〛 ≥ x.〚g〛`",
    /// where `g` is chosen such that `〚f*g〛` is a vector of `outbasis`.
    pub fn multiply_and_unlabel(self, outbasis: Basis<F>) -> Self
    where
        N: Neg<Output = N>,
    {
        assert_eq!(outbasis.t, Type::empty());
        let unlabeling = Unlabeling::total(self.meta.basis.t);
        let other = outbasis.with_type(self.meta.basis.t) / self.meta.basis;
        let splitcount = SplitCount::from_input(&self.meta.basis, &other);
        let operator = MulAndUnlabel {
            split: splitcount,
            unlabeling,
        };
        //
        let table: Vec<CsMat<N>> = operator
            .get()
            .iter()
            .map(|m| m.map(|&x| N::from_i64(x).unwrap()))
            .collect();
        //
        let mut data = Vec::new();
        //
        for ineq in self.data {
            ineq.multiply_by_all(&table, &mut data)
        }
        Self {
            data,
            meta: self.meta.multiply(other, Expr::Var(0)).untype(),
        }
    }
}

/// Return the vector corresponding to the unlabeled flag `f`.
pub fn flag<N, F>(f: &F) -> QFlag<N, F>
where
    N: Num + Clone,
    F: Flag,
{
    Basis::new(f.size()).flag(f)
}

/// Return the vector corresponding to the flag `f` with `type_size` labeled
/// vertices.
pub fn flag_typed<N, F>(f: &F, type_size: usize) -> QFlag<N, F>
where
    N: Num + Clone,
    F: Flag,
{
    let flag = f.canonical_typed(type_size);
    let type_flag = flag.induce(&(0..type_size).collect::<Vec<_>>()); // type
    let t = Type::from_flag(&type_flag);
    let basis = Basis::new(f.size()).with_type(t);
    basis.flag(&flag)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::Graph;
    use ndarray::array;
    #[test]
    fn test_internals() {
        // matching_scales
        assert_eq!(matching_scales(15, 12), (4, 5, 60));
        assert_eq!(matching_scales(2, 24), (12, 1, 24));
        let (c1, c2, scale): (u64, u64, _) = matching_scales(1788, 2444);
        let big = 1788 * 2444 * 1048;
        assert_eq!((big * c1) / scale, big / 1788);
        assert_eq!((big * c2) / scale, big / 2444);
    }
    #[test]
    fn test_qflags() {
        let qflag = QFlag {
            basis: Basis::<Graph>::new(1),
            data: array![3., 2., -5., 45.14],
            scale: 42,
            expr: Expr::Zero,
        };
        assert_eq!(qflag.clone().no_scale(), qflag)
    }
    #[test]
    fn test_qflag_pows() {
        let b = Basis::new(2).with_type(Type::new(1, 0));
        let v: QFlag<i128, Graph> = b.random();
        assert_eq!(v.pow(0), b.with_size(1).one());
        assert_eq!(v.pow(1), v);
        assert_eq!(v.pow(2), &v * &v);
        assert_eq!(v.pow(3), &v * &(&v * &v));
    }
}
