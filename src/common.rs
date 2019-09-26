//! Flat data structures for binary relations.

use crate::flags::Arc;
use crate::iterators;
use crate::iterators::StreamingIterator;
use std::fmt::Debug;
use std::iter::Chain;
use std::mem::swap;
use std::ops::{Index, IndexMut, Neg, Range};

/// Common interface for square matrices stored in a single array while
/// taking advantage of symetries.
pub trait FlatMatrix: Sized {
    /// Type of the entries of the matrix.
    type Item;
    /// Iterator on one line of the matrix.
    type LineIter: Iterator<Item = usize>;
    type IndexIter: Iterator<Item = (usize, usize)>;
    // needed
    /// Length of the underlying vector depending on the number of
    /// line of the matrix.
    fn data_size(size: usize) -> usize;
    /// Direct access to the underlying vector
    fn get_vec(&self) -> &Vec<Self::Item>;
    /// Mutable access to the underlying vector
    fn get_vec_mut(&mut self) -> &mut Vec<Self::Item>;
    /// Create a matrix by wrapping the underlying vector.
    fn from_vec(_: Vec<Self::Item>) -> Self;
    /// Map a pair of indices to the corresponding index in the matrix.
    fn flat_index(i: usize, j: usize) -> usize;
    /// Iterator on every non-symetric entries of the matrix.
    fn index_iter(n: usize) -> Self::IndexIter;
    /// Iterator on every non-symetric index an line `v`.
    fn line_iter(n: usize, v: usize) -> Self::LineIter;
    // recommended
    /// Give the number of line of the matrix
    fn possible_size(&self) -> usize {
        let mut res = 0;
        let data_size = self.get_vec().len();
        loop {
            debug_assert!(Self::data_size(res) <= data_size);
            if Self::data_size(res) == data_size {
                return res;
            } else {
                res += 1
            }
        }
    }
    // provided
    /// Access to a matrix element.
    #[inline]
    fn get(&self, i: usize, j: usize) -> Self::Item
    where
        Self::Item: Clone,
    {
        self.get_vec()[Self::flat_index(i, j)].clone()
    }
    /// Redefine a matrix entry.
    #[inline]
    fn set(&mut self, ij: (usize, usize), v: Self::Item) {
        self.get_vec_mut()[Self::flat_index(ij.0, ij.1)] = v
    }
    #[inline]
    /// Create a new matrix with `n` lines filled with `elem`.
    fn new(elem: Self::Item, n: usize) -> Self
    where
        Self::Item: Clone,
    {
        Self::from_vec(vec![elem; Self::data_size(n)])
    }
    /// Change the number of line in the matrix to `n`.
    /// If line are added, fill them with `e`.
    #[inline]
    fn resize(&mut self, n: usize, e: Self::Item)
    where
        Self::Item: Clone,
    {
        self.get_vec_mut().resize(Self::data_size(n), e)
    }
}

/// Relation R such that R(x,y) iff R(y,x).
#[derive(Clone, Debug, PartialOrd, Ord, Eq, PartialEq, Serialize, Deserialize)]
pub struct Sym<A>(Vec<A>);

use std;

impl<A> FlatMatrix for Sym<A> {
    type Item = A;
    type LineIter = Range<usize>;
    type IndexIter = Box<dyn Iterator<Item = (usize, usize)>>;
    #[inline]
    fn data_size(size: usize) -> usize {
        (size * (size + 1)) / 2
    }

    #[inline]
    fn flat_index(mut i: usize, mut j: usize) -> usize {
        if j < i {
            swap(&mut i, &mut j)
        };
        Self::data_size(j) + i
    }
    #[inline]
    fn get_vec(&self) -> &Vec<A> {
        &self.0
    }
    #[inline]
    fn get_vec_mut(&mut self) -> &mut Vec<A> {
        &mut self.0
    }
    #[inline]
    fn from_vec(v: Vec<A>) -> Self {
        Sym(v)
    }
    #[inline]
    fn index_iter(n: usize) -> Self::IndexIter {
        Box::new((0..n).flat_map(move |j| (0..=j).map(move |i| (i, j))))
    }
    #[inline]
    fn line_iter(n: usize, v: usize) -> Self::LineIter {
        debug_assert!(v <= n);
        (0..n)
    }
}

impl<A> Index<(usize, usize)> for Sym<A> {
    type Output = A;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.0[Self::flat_index(i, j)]
    }
}

impl<A> IndexMut<(usize, usize)> for Sym<A> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut A {
        &mut self.0[Self::flat_index(i, j)]
    }
}

/// Symetric relation R such that R(x,x) never holds
#[derive(Clone, Debug, PartialOrd, Ord, Eq, PartialEq, Serialize, Deserialize)]
pub struct SymNonRefl<A>(Vec<A>);

impl<A> FlatMatrix for SymNonRefl<A> {
    type Item = A;
    type LineIter = Chain<Range<usize>, Range<usize>>;
    type IndexIter = Box<dyn Iterator<Item = (usize, usize)>>;

    fn data_size(size: usize) -> usize {
        if size == 0 {
            0
        } else {
            ((size - 1) * (size)) / 2
        }
    }

    fn flat_index(mut i: usize, mut j: usize) -> usize {
        if j < i {
            swap(&mut i, &mut j)
        };
        debug_assert!(j > i);
        Self::data_size(j) + i
    }
    fn get_vec(&self) -> &Vec<A> {
        &self.0
    }
    fn get_vec_mut(&mut self) -> &mut Vec<A> {
        &mut self.0
    }
    fn from_vec(v: Vec<A>) -> Self {
        SymNonRefl(v)
    }
    fn index_iter(n: usize) -> Self::IndexIter {
        Box::new((0..n).flat_map(move |j| (0..j).map(move |i| (i, j))))
    }
    fn line_iter(n: usize, v: usize) -> Self::LineIter {
        debug_assert!(v <= n);
        (0..v).chain(v + 1..n)
    }
}

impl<A> Index<(usize, usize)> for SymNonRefl<A> {
    type Output = A;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.0[Self::flat_index(i, j)]
    }
}

impl<A> IndexMut<(usize, usize)> for SymNonRefl<A> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut A {
        &mut self.0[Self::flat_index(i, j)]
    }
}

/// Relation R such that R(x,y) = -R(y,x) and R(x,x) does not hold.
#[derive(Clone, Debug, PartialOrd, Ord, Eq, PartialEq, Serialize, Deserialize)]
pub struct AntiSym<A>(Vec<A>);

impl<A> AntiSym<A>
where
    A: Neg<Output = A> + Copy,
{
    fn flat_index_raw(i: usize, j: usize) -> usize {
        debug_assert!(j > i);
        Self::data_size(j) + i
    }
}

impl<A> FlatMatrix for AntiSym<A>
where
    A: Neg<Output = A> + Copy,
{
    type Item = A;
    type LineIter = Chain<Range<usize>, Range<usize>>;
    type IndexIter = Box<dyn Iterator<Item = (usize, usize)>>;

    fn data_size(size: usize) -> usize {
        SymNonRefl::<A>::data_size(size)
    }
    fn flat_index(i: usize, j: usize) -> usize {
        if j > i {
            Self::flat_index_raw(i, j)
        } else {
            Self::flat_index_raw(j, i)
        }
    }
    fn get_vec(&self) -> &Vec<A> {
        &self.0
    }
    fn get_vec_mut(&mut self) -> &mut Vec<A> {
        &mut self.0
    }
    fn from_vec(v: Vec<A>) -> Self {
        AntiSym(v)
    }
    fn index_iter(n: usize) -> Self::IndexIter {
        Box::new((0..n).flat_map(move |i| (0..i).map(move |j| (i, j))))
    }
    fn line_iter(n: usize, v: usize) -> Self::LineIter {
        debug_assert!(v <= n);
        (0..v).chain(v + 1..n)
    }
    fn get(&self, i: usize, j: usize) -> A {
        debug_assert!(i != j);
        if i < j {
            self.0[Self::flat_index_raw(i, j)]
        } else {
            -self.0[Self::flat_index_raw(j, i)]
        }
    }
    fn set(&mut self, (i, j): (usize, usize), v: A) {
        if i < j {
            self.0[Self::flat_index_raw(i, j)] = v
        } else {
            self.0[Self::flat_index_raw(j, i)] = -v
        }
    }
}

/// A trait that gives access to the list of the possible values
/// of a type.
pub trait Enum: Sized {
    /// List of possible values of `Self`.
    fn variants() -> Vec<Self>;
}

impl Enum for bool {
    fn variants() -> Vec<Self> {
        vec![true, false]
    }
}

impl Enum for Arc {
    fn variants() -> Vec<Self> {
        vec![Arc::Edge, Arc::BackEdge, Arc::None]
    }
}

/// Generic trait for binary relations.
pub trait BinRelation: Ord + Debug + Clone {
    fn invariant(&self, v: usize) -> Vec<Vec<usize>>;
    fn invariant_size() -> usize;
    fn induce(&self, p: &[usize]) -> Self;
    fn empty() -> Self;
    fn extensions(&self, n: usize) -> Vec<Self>;
}

impl<S> BinRelation for S
where
    S: FlatMatrix + Debug + Clone + Ord,
    S::Item: Enum + Ord + Clone + Copy + Debug + Default,
{
    fn invariant(&self, v: usize) -> Vec<Vec<usize>> {
        let mut variants = S::Item::variants();
        variants.sort();
        let mut res: Vec<Vec<usize>> = vec![Vec::new(); variants.len() - 1];
        let n = self.possible_size();
        for u in Self::line_iter(n, v) {
            let val = self.get(u, v);
            let e = variants.iter().position(|x| x == &val).unwrap();
            if e < variants.len() - 1 {
                res[e].push(u);
            }
        }
        res
    }
    fn invariant_size() -> usize {
        S::Item::variants().len() - 1
    }
    fn extensions(&self, n: usize) -> Vec<Self> {
        assert_eq!(self.get_vec().len(), Self::data_size(n));
        let mut res = Vec::new();
        let variants = S::Item::variants();
        let line_size = Self::line_iter(n + 1, n).count();
        let mut iter = iterators::Functions::new(line_size, variants.len());
        while let Some(f) = iter.next() {
            let mut new = self.clone();
            new.resize(n + 1, S::Item::default());
            for v in Self::line_iter(n + 1, n) {
                new.get_vec_mut()[Self::flat_index(v, n)] = variants[f[v]];
            }
            res.push(new);
        }
        res
    }
    fn induce(&self, p: &[usize]) -> Self
    where
        S::Item: Default + Clone,
    {
        let n = p.len();
        let mut res = Self::new(S::Item::default(), n);
        for (u1, u2) in Self::index_iter(n) {
            res.set((u1, u2), self.get(p[u1], p[u2]));
        }
        res
    }
    fn empty() -> Self {
        Self::new(S::Item::default(), 0)
    }
}

// Extension of BinRelation to small tuples: should be automatized
// size 2
impl<R1, R2> BinRelation for (R1, R2)
where
    R1: BinRelation,
    R2: BinRelation,
{
    fn invariant(&self, v: usize) -> Vec<Vec<usize>> {
        let mut res = self.0.invariant(v);
        res.append(&mut self.1.invariant(v));
        res
    }
    fn invariant_size() -> usize {
        R1::invariant_size() + R2::invariant_size()
    }
    fn induce(&self, p: &[usize]) -> Self {
        (self.0.induce(p), self.1.induce(p))
    }
    fn empty() -> Self {
        (R1::empty(), R2::empty())
    }
    fn extensions(&self, n: usize) -> Vec<Self> {
        let e1 = self.0.extensions(n);
        let e2 = self.1.extensions(n);
        let mut res = Vec::new();
        for x1 in &e1 {
            for x2 in &e2 {
                res.push((x1.clone(), x2.clone()))
            }
        }
        res
    }
}
// size 3
impl<R1, R2, R3> BinRelation for (R1, R2, R3)
where
    R1: BinRelation,
    R2: BinRelation,
    R3: BinRelation,
{
    fn invariant(&self, v: usize) -> Vec<Vec<usize>> {
        let mut res = self.0.invariant(v);
        res.append(&mut self.1.invariant(v));
        res.append(&mut self.2.invariant(v));
        res
    }
    fn invariant_size() -> usize {
        R1::invariant_size() + R2::invariant_size() + R3::invariant_size()
    }
    fn induce(&self, p: &[usize]) -> Self {
        (self.0.induce(p), self.1.induce(p), self.2.induce(p))
    }
    fn empty() -> Self {
        (R1::empty(), R2::empty(), R3::empty())
    }
    fn extensions(&self, n: usize) -> Vec<Self> {
        let e1 = self.0.extensions(n);
        let e2 = self.1.extensions(n);
        let e3 = self.2.extensions(n);
        let mut res = Vec::new();
        for x1 in &e1 {
            for x2 in &e2 {
                for x3 in &e3 {
                    res.push((x1.clone(), x2.clone(), x3.clone()))
                }
            }
        }
        res
    }
}

/// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::*;
    use crate::Flag;

    fn auto_test_flatmatrix<S: FlatMatrix<Item = i64>>(n: usize) {
        assert_eq!(S::index_iter(n).count(), S::data_size(n));
        println!("{}_{}", n + 1, n / 2);
        assert_eq!(
            S::line_iter(n + 1, n / 2).count(),
            S::data_size(n + 1) - S::data_size(n)
        );
        // injectivity of line_iter
        for u in 0..n {
            let mut x = S::new(0, n);
            for v in S::line_iter(n, u) {
                assert_eq!(x.get(u, v), 0);
                x.set((u, v), 1);
            }
        }
        // injectivity of index_iter
        let mut x = S::new(0, n);
        for (u, v) in S::index_iter(n) {
            assert_eq!(x.get(u, v), 0);
            x.set((u, v), 1);
        }
    }

    #[test]
    fn generate_graph() {
        for (size, &nb) in [1, 1, 2, 4, 11, 34].iter().enumerate() {
            assert_eq!(Graph::generate(size).len(), nb);
        }
    }

    // #[test]
    // fn symflatmatrix_size_unit() {
    //     for i in 0..4 {
    //         let m = Sym::new(0,i);
    //         assert_eq!(i, m.size())
    //     }
    // }

    #[test]
    fn symnonrefl() {
        assert_eq!(SymNonRefl::new(42, 0).0.len(), 0);
        assert_eq!(SymNonRefl::new(42, 5)[(4, 3)], 42);
        let mut m = SymNonRefl::new(0, 12);
        m[(3, 2)] = 11;
        assert_eq!(m[(2, 3)], 11);
        m[(3, 4)] = 22;
        assert_eq!(m[(4, 3)], 22);

        let n = 10;
        let mut m = SymNonRefl::new(0, n);
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    m[(i, j)] += 1;
                }
            }
        }
        for &x in m.0.iter() {
            assert_eq!(x, 2)
        }
        //        assert_eq!(SymNonRefl::<bool>::invariant_size(), 1)
    }

    #[test]
    fn antisym() {
        assert_eq!(AntiSym::new(42, 0).0.len(), 0);
        let mut rel = AntiSym::new(0, 12);
        rel.set((5, 2), 42);
        assert_eq!(rel.get(5, 2), 42);
        assert_eq!(rel.get(2, 5), -42);
        let n = 10;
        let mut m = AntiSym::new(0, n);
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let v = m.get(i, j);
                    if i < j {
                        assert_eq!(v, 0);
                        m.set((i, j), 42)
                    } else {
                        assert_eq!(v, -42)
                    }
                }
            }
        }
        for &x in m.0.iter() {
            assert!(x != 0)
        }
    }

    #[test]
    fn flatmatrix_generic() {
        auto_test_flatmatrix::<Sym<_>>(0);
        auto_test_flatmatrix::<Sym<_>>(1);
        auto_test_flatmatrix::<Sym<_>>(5);
        //
        auto_test_flatmatrix::<AntiSym<_>>(0);
        auto_test_flatmatrix::<AntiSym<_>>(1);
        auto_test_flatmatrix::<AntiSym<_>>(5);
        //
        // auto_test_flatmatrix::<NonRefl<_>>(0);
        // auto_test_flatmatrix::<NonRefl<_>>(1);
        // auto_test_flatmatrix::<NonRefl<_>>(5);
        //
        auto_test_flatmatrix::<SymNonRefl<_>>(0);
        auto_test_flatmatrix::<SymNonRefl<_>>(1);
        auto_test_flatmatrix::<SymNonRefl<_>>(5);
    }
}
