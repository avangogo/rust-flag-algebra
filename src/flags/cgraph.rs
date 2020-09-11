//! Example of flags: Graphs.

use crate::combinatorics;
use crate::common::*;
use crate::flag::{Flag, SubClass};
use crate::flags::Colored;
use crate::iterators;
use crate::iterators::StreamingIterator;
use canonical_form::Canonize;
use core::marker::PhantomData;
use std::fmt;
use std::fmt::Debug;
use typenum::Unsigned;

/// An undirected graph.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
pub struct CGraph<E> {
    pub size: usize,
    edge: SymNonRefl<u8>,
    phantom: PhantomData<E>,
}

impl<E> CGraph<E> {
    /// Returns `true` id `uv` is an edge.
    pub fn is_edge(&self, u: usize, v: usize) -> bool {
        u != v && self.edge[(u, v)] > 0
    }

    pub fn edge(&self, u: usize, v: usize) -> u8 {
        self.edge[(u, v)]
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

    pub fn new(n: usize, edge: &[((usize, usize), u8)]) -> Self
    where
        E: Unsigned,
    {
        let mut res = Self::empty(n);
        for &((u, v), val) in edge {
            assert!(u < n);
            assert!(v < n);
            assert!(u != v);
            assert!(val < E::U8);
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
            phantom: PhantomData,
        }
    }
}

impl<E> fmt::Display for CGraph<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let sep = if self.size < 10 { "" } else { "-" };
        write!(f, "(|V|={}, E={{", self.size)?;
        for u in 0..self.size {
            for v in 0..u {
                let uv = self.edge[(u, v)];
                if uv > 0 {
                    write!(f, " {}{}{}:{}", v, sep, u, uv)?
                }
            }
        }
        write!(f, " }})")
    }
}

impl<E> Canonize for CGraph<E>
where
    E: Unsigned + Clone + Eq + Ord + Debug,
{
    fn size(&self) -> usize {
        self.size
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        assert!(v < self.size);
        let mut res = vec![Vec::new(); E::USIZE - 1];
        for u in 0..self.size() {
            if u != v {
                let edge = self.edge[(u, v)];
                if edge > 0 {
                    res[edge as usize - 1].push(u)
                }
            }
        }
        res
    }

    fn apply_morphism(&self, p: &[usize]) -> Self {
        self.induce(&combinatorics::invert(&p))
    }
}

impl<E> Flag for CGraph<E>
where
    E: Unsigned + Clone + Eq + Ord + Debug,
{
    fn induce(&self, p: &[usize]) -> Self {
        Self {
            size: p.len(),
            edge: self.edge.induce0(p),
            phantom: PhantomData,
        }
    }

    const NAME: &'static str = "CGraph_FIXME";

    fn size_zero_flags() -> Vec<Self> {
        vec![Self::empty(0)]
    }

    fn superflags(&self) -> Vec<Self> {
        let n = self.size;
        let mut res = Vec::new();
        let mut iter = iterators::Functions::new(n, E::USIZE);
        while let Some(map) = iter.next() {
            let mut edge = self.edge.clone();
            edge.resize(n + 1, 0);
            for (v, &c) in map.iter().enumerate() {
                edge[(v, n)] = c as u8;
            }
            res.push(Self {
                edge,
                size: n + 1,
                phantom: PhantomData,
            });
        }
        res
    }
}

/// This should be properly implemented with a trait:

impl<A, E> SubClass<CGraph<E>, A> {
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

impl<A, E> Colored<CGraph<E>, A> {
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

impl<A, B, E> SubClass<Colored<CGraph<E>, A>, B> {
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
    use typenum::{U1, U2, U3};

    #[test]
    fn new_unit() {
        let g =
            CGraph::<U3>::new(4, &[((0, 1), 1), ((2, 3), 1), ((0, 2), 2), ((0, 3), 2)]).canonical();
        let h =
            CGraph::<U3>::new(4, &[((1, 0), 2), ((2, 1), 2), ((3, 1), 1), ((0, 2), 1)]).canonical();
        assert_eq!(g, h);
    }

    #[test]
    fn generate_cgraph() {
        type G1 = CGraph<U1>;
        type G2 = CGraph<U2>;
        type G3 = CGraph<U3>;

        assert_eq!(G1::generate(12).len(), 1);
        assert_eq!(G2::generate(6).len(), 156);
        assert_eq!(G3::generate(3).len(), 10);
    }
}
