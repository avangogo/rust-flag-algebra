//! Example of flags: Graphs.

use crate::flag::{Flag, SubClass};
use crate::flags::Colored;
use crate::flags::common::*;
use crate::iterators;
use crate::iterators::StreamingIterator;
use canonical_form::Canonize;
use std::fmt;
use std::fmt::Debug;

/// Edge-colored graphs.
///
/// Implentation of graphs whose edges are colored by a number in a set `{1,…,K-1}`
/// (So that accounting for non-edges, each pair has `K` possible states).
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
pub struct CGraph<const K: u8> {
    pub size: usize,
    edge: SymNonRefl<u8>,
}

impl<const K: u8> CGraph<K> {
    /// Returns `true` if `uv` is an edge.
    pub fn is_edge(&self, u: usize, v: usize) -> bool {
        self.edge.get(u, v).is_some()
    }
    /// Return the color of an edge `uv`.
    /// Returns `0` if there is no edge, `u == v` included.
    pub fn edge(&self, u: usize, v: usize) -> u8 {
        self.edge.get(u, v).unwrap_or(0)
    }
    /// Return the vector of vertices adjacent to `v`.
    pub fn nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if self.is_edge(u, v) {
                res.push(u);
            }
        }
        res
    }
    /// Create a new colored graph with the specified edges.
    /// Edges not listed get value 0.
    pub fn new(n: usize, edge: &[((usize, usize), u8)]) -> Self {
        let mut res = Self::empty(n);
        for &((u, v), val) in edge {
            assert!(u < n);
            assert!(v < n);
            assert!(u != v);
            assert!(val < K);
            assert_eq!(res.edge[(u, v)], 0);
            res.edge[(u, v)] = val;
        }
        res
    }
    /// Create the graph on `n` vertices with no edge.  
    pub fn empty(n: usize) -> Self {
        Self {
            size: n,
            edge: SymNonRefl::new(0, n),
        }
    }
}

impl<const K: u8> fmt::Display for CGraph<K> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let sep = if self.size < 10 { "" } else { "-" };
        write!(f, "(|V|={}, E={{", self.size)?;
        for u in 0..self.size {
            for v in 0..u {
                let uv = self.edge[(u, v)];
                if uv > 0 {
                    write!(f, " {v}{sep}{u}:{uv}")?
                }
            }
        }
        write!(f, " }})")
    }
}

impl<const K: u8> Canonize for CGraph<K> {
    fn size(&self) -> usize {
        self.size
    }
    fn invariant_neighborhood(&self, v: usize) -> impl Iterator<Item = (usize, u64)> {
        assert!(v < self.size);
        // `cell`, not `get`: the hottest read in the library, and `u != v` is
        // already tested here.
        (0..self.size())
            .filter(move |&u| u != v)
            .map(move |u| (u, self.edge.cell(u, v) as u64))
            .filter(|(_, weight)| *weight > 0)
    }

    fn apply_morphism(&self, p: &[usize]) -> Self {
        Self {
            size: self.size,
            edge: self.edge.relabelled(p),
        }
    }
}

impl<const K: u8> Flag for CGraph<K> {
    fn induce(&self, p: &[usize]) -> Self {
        Self {
            size: p.len(),
            edge: self.edge.induced(p),
        }
    }

    fn name() -> String {
        format!("CGraph_{K}")
    }

    fn size_zero_flags() -> Vec<Self> {
        vec![Self::empty(0)]
    }

    fn superflags(&self) -> Vec<Self> {
        let n = self.size;
        let mut res = Vec::new();
        let mut iter = iterators::Functions::new(n, K as usize);
        while let Some(map) = iter.next() {
            let mut edge = self.edge.clone();
            edge.resize(n + 1, 0);
            for (v, &c) in map.iter().enumerate() {
                edge[(v, n)] = c as u8;
            }
            res.push(Self { edge, size: n + 1 });
        }
        res
    }
}

/// This should be properly implemented with a trait:
impl<A, const K: u8> SubClass<CGraph<K>, A> {
    pub fn size(&self) -> usize {
        self.content.size
    }
    pub fn is_edge(&self, u: usize, v: usize) -> bool {
        self.content.is_edge(u, v)
    }
    pub fn edge(&self, u: usize, v: usize) -> u8 {
        self.content.edge(u, v)
    }
}

impl<const K: u8, const N: u8> Colored<CGraph<K>, N> {
    pub fn size(&self) -> usize {
        self.content.size
    }
    pub fn is_edge(&self, u: usize, v: usize) -> bool {
        self.content.is_edge(u, v)
    }
    pub fn edge(&self, u: usize, v: usize) -> u8 {
        self.content.edge(u, v)
    }
}

impl<A, const K: u8, const N: u8> SubClass<Colored<CGraph<K>, N>, A> {
    pub fn size(&self) -> usize {
        self.content.size()
    }
    pub fn is_edge(&self, u: usize, v: usize) -> bool {
        self.content.is_edge(u, v)
    }
    pub fn edge(&self, u: usize, v: usize) -> u8 {
        self.content.edge(u, v)
    }
}

/// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_is_not_capped() {
        assert_eq!(<CGraph<3> as Flag>::name(), "CGraph_3");
        // The old macro array only had entries for K <= 8.
        assert_eq!(<CGraph<12> as Flag>::name(), "CGraph_12");
    }

    /// `edge(v, v)` used to read the first cell of the next row: here it
    /// answered 2 for `edge(1, 1)`, and `edge(3, 3)` panicked.
    #[test]
    fn the_diagonal_is_not_an_edge() {
        let g = CGraph::<3>::new(4, &[((0, 1), 1), ((0, 2), 2)]);
        for v in 0..4 {
            assert_eq!(g.edge(v, v), 0);
            assert!(!g.is_edge(v, v));
        }
    }

    #[test]
    fn new_unit() {
        let g =
            CGraph::<3>::new(4, &[((0, 1), 1), ((2, 3), 1), ((0, 2), 2), ((0, 3), 2)]).canonical();
        let h =
            CGraph::<3>::new(4, &[((1, 0), 2), ((2, 1), 2), ((3, 1), 1), ((0, 2), 1)]).canonical();
        assert_eq!(g, h);
    }

    #[test]
    fn generate_cgraph() {
        type G1 = CGraph<1>;
        type G2 = CGraph<2>;
        type G3 = CGraph<3>;

        assert_eq!(G1::generate(12).len(), 1);
        assert_eq!(G2::generate(6).len(), 156);
        assert_eq!(G3::generate(3).len(), 10);
    }
}
