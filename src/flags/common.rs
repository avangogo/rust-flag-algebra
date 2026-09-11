//! Flat data structures for binary relations.

use crate::flags::Arc;
use std::fmt::Debug;
use std::mem::swap;
use std::ops::{Index, IndexMut, Neg, Range};

/// Common interface for square matrices stored in a single array while
/// taking advantage of symetries.
pub trait FlatMatrix: Sized {
    /// Type of the entries of the matrix.
    type Item;
    // needed
    /// Length of the underlying vector depending on the number of
    /// line of the matrix.
    fn data_size(size: usize) -> usize;
    /// Direct access to the underlying vector
    fn data(&self) -> &[Self::Item];
    /// Mutable access to the underlying vector
    fn data_mut(&mut self) -> &mut Vec<Self::Item>;
    /// Create a matrix by wrapping the underlying vector.
    fn from_vec(_: Vec<Self::Item>) -> Self;
    /// Map a pair of indices to the corresponding index in the matrix.
    fn flat_index(i: usize, j: usize) -> usize;
    /// Iterator on every non-symetric index an line `v`.
    fn halfline_iter(v: usize) -> Range<usize>;
    // provided
    /// Access to a matrix element.
    #[inline]
    fn cell(&self, i: usize, j: usize) -> Self::Item
    where
        Self::Item: Clone,
    {
        self.data()[Self::flat_index(i, j)].clone()
    }
    /// Redefine a matrix entry.
    #[inline]
    fn set_cell(&mut self, ij: (usize, usize), v: Self::Item) {
        self.data_mut()[Self::flat_index(ij.0, ij.1)] = v
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
        self.data_mut().resize(Self::data_size(n), e)
    }
}

/// The cells a flat matrix holds, as a stream.
///
/// Each implementation fixes its own protocol: which ordered pairs name a cell,
/// and what one carries. `AntiSym<Arc>` names a cell by an ordered pair, so a
/// reciprocal pair is two cells that `build` folds back into one.
///
/// `build(n, m.iter())` returns `m`. The order is the caller's to supply:
/// a matrix of order 0 and one of order 1 both hold no cell.
pub trait Relation: Sized {
    /// What a cell carries; `()` when being there is all there is to it.
    type Value;
    /// The ordered pairs that name a cell of a matrix of order `n`, each once.
    fn positions(n: usize) -> impl Iterator<Item = (usize, usize)>;
    /// The content of cell `(u, v)`, or `None` when it is not there.
    ///
    /// Defined on every ordered pair, `u == v` included; under a symmetric
    /// protocol that is more pairs than [`Relation::iter`] yields.
    fn get(&self, u: usize, v: usize) -> Option<Self::Value>;
    /// Every cell that is there, once, in this implementation's protocol.
    fn iter(&self) -> impl Iterator<Item = (usize, usize, Self::Value)>;
    /// The matrix of order `n` holding `entries`; every cell the stream omits
    /// is absent.
    fn build(n: usize, entries: impl IntoIterator<Item = (usize, usize, Self::Value)>) -> Self;

    // provided
    /// The submatrix on `p`, where vertex `i` of the result is `p[i]`.
    fn induced(&self, p: &[usize]) -> Self {
        Self::build(
            p.len(),
            Self::positions(p.len())
                .filter_map(|(u, v)| self.get(p[u], p[v]).map(|value| (u, v, value))),
        )
    }
    /// `self` with vertex `i` renamed to `p[i]`, for a permutation `p` of the
    /// vertices.
    fn relabelled(&self, p: &[usize]) -> Self {
        Self::build(
            p.len(),
            self.iter().map(|(u, v, value)| (p[u], p[v], value)),
        )
    }
}

/// The pairs `(u, v)` with `u < v < n`.
fn pairs(n: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..n).flat_map(move |v| (0..v).map(move |u| (u, v)))
}

/// The pairs `(u, v)`, `u < v`, in the order `SymNonRefl` and `AntiSym` store
/// them, so that walking the data needs no `flat_index`.
fn stored_pairs(len: usize) -> impl Iterator<Item = (usize, usize)> {
    let (mut u, mut v) = (0, 1);
    (0..len).map(move |_| {
        let pair = (u, v);
        u += 1;
        if u == v {
            u = 0;
            v += 1;
        }
        pair
    })
}

/// The pairs `(u, v)` with `u != v` and both below `n`.
fn ordered_pairs(n: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..n).flat_map(move |v| (0..n).filter(move |&u| u != v).map(move |u| (u, v)))
}

impl Relation for SymNonRefl<bool> {
    type Value = ();
    fn positions(n: usize) -> impl Iterator<Item = (usize, usize)> {
        pairs(n)
    }
    #[inline]
    fn get(&self, u: usize, v: usize) -> Option<()> {
        (u != v && self.cell(u, v)).then_some(())
    }
    fn iter(&self) -> impl Iterator<Item = (usize, usize, ())> {
        stored_pairs(self.data().len())
            .zip(self.data())
            .filter_map(|((u, v), &present)| present.then_some((u, v, ())))
    }
    fn build(n: usize, entries: impl IntoIterator<Item = (usize, usize, ())>) -> Self {
        let mut res = Self::new(false, n);
        for (u, v, ()) in entries {
            debug_assert!(u != v && u < n && v < n);
            res.set_cell((u, v), true);
        }
        res
    }
    /// Moves absent cells too: the write is a no-op, and the test that would
    /// skip it is one the branch predictor cannot call.
    fn relabelled(&self, p: &[usize]) -> Self {
        let mut res = Self::new(false, p.len());
        for ((u, v), &present) in stored_pairs(self.data().len()).zip(self.data()) {
            res.set_cell((p[u], p[v]), present);
        }
        res
    }
}

impl Relation for SymNonRefl<u8> {
    type Value = u8;
    fn positions(n: usize) -> impl Iterator<Item = (usize, usize)> {
        pairs(n)
    }
    #[inline]
    fn get(&self, u: usize, v: usize) -> Option<u8> {
        if u == v {
            None
        } else {
            match self.cell(u, v) {
                0 => None,
                color => Some(color),
            }
        }
    }
    fn iter(&self) -> impl Iterator<Item = (usize, usize, u8)> {
        stored_pairs(self.data().len())
            .zip(self.data())
            .filter_map(|((u, v), &color)| (color > 0).then_some((u, v, color)))
    }
    fn build(n: usize, entries: impl IntoIterator<Item = (usize, usize, u8)>) -> Self {
        let mut res = Self::new(0, n);
        for (u, v, color) in entries {
            debug_assert!(u != v && u < n && v < n);
            debug_assert!(color > 0, "0 is absence, not a colour");
            res.set_cell((u, v), color);
        }
        res
    }
    fn relabelled(&self, p: &[usize]) -> Self {
        let mut res = Self::new(0, p.len());
        for ((u, v), &color) in stored_pairs(self.data().len()).zip(self.data()) {
            res.set_cell((p[u], p[v]), color);
        }
        res
    }
}

impl Relation for AntiSym<Arc> {
    type Value = ();
    fn positions(n: usize) -> impl Iterator<Item = (usize, usize)> {
        ordered_pairs(n)
    }
    #[inline]
    fn get(&self, u: usize, v: usize) -> Option<()> {
        if u == v {
            None
        } else {
            matches!(self.cell(u, v), Arc::Edge | Arc::Reciprocal).then_some(())
        }
    }
    fn iter(&self) -> impl Iterator<Item = (usize, usize, ())> {
        stored_pairs(self.data().len())
            .zip(self.data())
            .flat_map(|((u, v), &arc)| {
                let (forward, backward) = match arc {
                    Arc::Edge => (true, false),
                    Arc::BackEdge => (false, true),
                    Arc::Reciprocal => (true, true),
                    Arc::None => (false, false),
                };
                [
                    forward.then_some((u, v, ())),
                    backward.then_some((v, u, ())),
                ]
            })
            .flatten()
    }
    fn build(n: usize, entries: impl IntoIterator<Item = (usize, usize, ())>) -> Self {
        let mut res = Self::new(Arc::None, n);
        for (u, v, ()) in entries {
            debug_assert!(u != v && u < n && v < n);
            // `cell` is oriented: this is the arc `u -> v`.
            let merged = match res.cell(u, v) {
                Arc::None => Arc::Edge,
                Arc::BackEdge => Arc::Reciprocal,
                already => already,
            };
            res.set_cell((u, v), merged);
        }
        res
    }
    /// One `Arc` answers for both directions, so the cells can be gathered
    /// whole. Order matters: reading the pair as `(p[j], p[i])` reverses it.
    fn induced(&self, p: &[usize]) -> Self {
        let mut res = Self::new(Arc::None, p.len());
        for (j, &pj) in p.iter().enumerate() {
            for (i, &pi) in p[..j].iter().enumerate() {
                res.set_cell((i, j), self.cell(pi, pj));
            }
        }
        res
    }
    /// Moving whole cells is half the writes of the generic fold, and no reads.
    fn relabelled(&self, p: &[usize]) -> Self {
        let mut res = Self::new(Arc::None, p.len());
        for ((u, v), &arc) in stored_pairs(self.data().len()).zip(self.data()) {
            // `set_cell` negates when `p[u] > p[v]`, and `-None` is `None`.
            res.set_cell((p[u], p[v]), arc);
        }
        res
    }
}

/// Relation R such that R(x,y) iff R(y,x).
#[derive(Clone, Debug, PartialOrd, Ord, Eq, PartialEq, Serialize, Deserialize)]
pub struct Sym<A>(Vec<A>);

impl<A> FlatMatrix for Sym<A> {
    type Item = A;
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
    fn data(&self) -> &[A] {
        &self.0
    }
    #[inline]
    fn data_mut(&mut self) -> &mut Vec<A> {
        &mut self.0
    }
    #[inline]
    fn from_vec(v: Vec<A>) -> Self {
        Sym(v)
    }
    #[inline]
    #[allow(clippy::range_plus_one)]
    fn halfline_iter(v: usize) -> Range<usize> {
        0..(v + 1)
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
    fn data(&self) -> &[A] {
        &self.0
    }
    fn data_mut(&mut self) -> &mut Vec<A> {
        &mut self.0
    }
    fn from_vec(v: Vec<A>) -> Self {
        SymNonRefl(v)
    }
    #[inline]
    fn halfline_iter(v: usize) -> Range<usize> {
        0..v
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
    fn data(&self) -> &[A] {
        &self.0
    }
    fn data_mut(&mut self) -> &mut Vec<A> {
        &mut self.0
    }
    fn from_vec(v: Vec<A>) -> Self {
        AntiSym(v)
    }
    fn cell(&self, i: usize, j: usize) -> A {
        debug_assert!(i != j);
        if i < j {
            self.0[Self::flat_index_raw(i, j)]
        } else {
            -self.0[Self::flat_index_raw(j, i)]
        }
    }
    fn set_cell(&mut self, (i, j): (usize, usize), v: A) {
        if i < j {
            self.0[Self::flat_index_raw(i, j)] = v
        } else {
            self.0[Self::flat_index_raw(j, i)] = -v
        }
    }
    #[inline]
    fn halfline_iter(v: usize) -> Range<usize> {
        0..v
    }
}

/// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Flag;
    use crate::combinatorics::invert;
    use crate::flags::*;
    use crate::iterators;
    use crate::iterators::StreamingIterator;

    /// `halfline_iter` names the cells of each line; `flat_index` must map
    /// them one to one onto the data, and ignore the order of the pair.
    fn auto_test_flatmatrix<S: FlatMatrix<Item = i64>>(n: usize) {
        let mut seen = vec![false; S::data_size(n)];
        for j in 0..n {
            for i in S::halfline_iter(j) {
                let k = S::flat_index(i, j);
                assert_eq!(S::flat_index(j, i), k, "flat_index reads the order");
                assert!(!seen[k], "two cells share the index {k}");
                seen[k] = true;
            }
        }
        assert!(seen.into_iter().all(|b| b), "some cell is unreachable");
    }

    #[test]
    fn generate_graph() {
        for (size, &nb) in [1, 1, 2, 4, 11, 34].iter().enumerate() {
            assert_eq!(Graph::generate(size).len(), nb);
        }
    }

    /// Every matrix of order `n` over `alphabet`.
    fn every_matrix<S: FlatMatrix>(n: usize, alphabet: &[S::Item]) -> Vec<S>
    where
        S::Item: Clone,
    {
        let cells = S::data_size(n);
        (0..alphabet.len().pow(cells as u32))
            .map(|mut code| {
                S::from_vec(
                    (0..cells)
                        .map(|_| {
                            let digit = code % alphabet.len();
                            code /= alphabet.len();
                            alphabet[digit].clone()
                        })
                        .collect(),
                )
            })
            .collect()
    }

    /// The gather each flag class hand-rolled before [`Relation`].
    fn gather<S: FlatMatrix>(m: &S, p: &[usize], empty: S::Item) -> S
    where
        S::Item: Clone,
    {
        let mut res = S::new(empty, p.len());
        for u1 in 0..p.len() {
            for u2 in 0..u1 {
                res.set_cell((u1, u2), m.cell(p[u1], p[u2]));
            }
        }
        res
    }

    /// Everything a protocol owes: the round trip, agreement between `get` and
    /// `iter`, and `induced`/`relabelled` against the gather they replace.
    fn check_protocol<S, V>(n: usize, m: &S, empty: S::Item)
    where
        S: Relation<Value = V> + FlatMatrix + PartialEq + Debug,
        S::Item: Clone,
        V: PartialEq + Debug + Copy,
    {
        let cells: Vec<_> = m.iter().collect();
        assert_eq!(S::build(n, cells.clone()), *m, "round trip");

        for u in 0..n {
            assert_eq!(m.get(u, u), None, "the diagonal is never a cell");
        }
        let mut from_iter = cells;
        let mut from_get: Vec<_> = S::positions(n)
            .filter_map(|(u, v)| m.get(u, v).map(|value| (u, v, value)))
            .collect();
        from_iter.sort_by_key(|&(u, v, _)| (u, v));
        from_get.sort_by_key(|&(u, v, _)| (u, v));
        assert_eq!(from_iter, from_get, "get and iter disagree");

        let mut perms = iterators::Injection::new(n, n);
        while let Some(p) = perms.next() {
            assert_eq!(m.induced(p), gather(m, p, empty.clone()), "induced");
            assert_eq!(
                m.relabelled(p),
                gather(m, &invert(p), empty.clone()),
                "relabelled"
            );
        }
    }

    #[test]
    fn protocols() {
        const ARCS: [Arc; 4] = [Arc::None, Arc::Edge, Arc::BackEdge, Arc::Reciprocal];
        for n in 0..5 {
            for m in every_matrix::<SymNonRefl<bool>>(n, &[false, true]) {
                check_protocol(n, &m, false);
            }
            for m in every_matrix::<SymNonRefl<u8>>(n, &[0, 1, 2]) {
                check_protocol(n, &m, 0);
            }
            for m in every_matrix::<AntiSym<Arc>>(n, &ARCS) {
                check_protocol(n, &m, Arc::None);
            }
        }
    }

    /// A reciprocal pair is two cells, and `build` folds them back into one.
    #[test]
    fn arcs_are_one_cell_per_direction() {
        let mut m = AntiSym::new(Arc::None, 3);
        m.set_cell((0, 1), Arc::Edge);
        m.set_cell((1, 2), Arc::Reciprocal);
        assert_eq!(
            m.iter().collect::<Vec<_>>(),
            vec![(0, 1, ()), (1, 2, ()), (2, 1, ())]
        );
        let fold = |entries: Vec<(usize, usize, ())>| AntiSym::<Arc>::build(3, entries).cell(1, 2);
        assert_eq!(fold(vec![(2, 1, ()), (1, 2, ())]), Arc::Reciprocal);
        assert_eq!(fold(vec![(1, 2, ()), (1, 2, ())]), Arc::Edge);
    }

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
        for &x in &m.0 {
            assert_eq!(x, 2)
        }
    }

    #[test]
    fn antisym() {
        assert_eq!(AntiSym::new(42, 0).0.len(), 0);
        let mut rel = AntiSym::new(0, 12);
        rel.set_cell((5, 2), 42);
        assert_eq!(rel.cell(5, 2), 42);
        assert_eq!(rel.cell(2, 5), -42);
        let n = 10;
        let mut m = AntiSym::new(0, n);
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let v = m.cell(i, j);
                    if i < j {
                        assert_eq!(v, 0);
                        m.set_cell((i, j), 42)
                    } else {
                        assert_eq!(v, -42)
                    }
                }
            }
        }
        for &x in &m.0 {
            assert!(x != 0)
        }
    }

    #[test]
    fn flatmatrix_generic() {
        for n in [0, 1, 5] {
            auto_test_flatmatrix::<Sym<_>>(n);
            auto_test_flatmatrix::<AntiSym<_>>(n);
            auto_test_flatmatrix::<SymNonRefl<_>>(n);
        }
    }
}
