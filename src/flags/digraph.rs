//! Example of flags: directed graphs.
use crate::combinatorics;
use crate::common::*;
use crate::flag::{Flag, SubClass, SubFlag};
use crate::iterators::{Functions, StreamingIterator};
use canonical_form::Canonize;
use std::fmt;
use std::ops::Neg;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy, Serialize, Deserialize)]
/// The arc between two nodes of a directed graphs.
pub enum Arc {
    /// No arc.
    None,
    /// Arc from the first to the second vertex considered.
    Edge,
    /// Arc from the second to the first vertex considered.
    BackEdge,
}

impl Default for Arc {
    fn default() -> Self {
        None
    }
}

impl Neg for Arc {
    type Output = Self;

    fn neg(self) -> Self {
        match self {
            Edge => BackEdge,
            BackEdge => Edge,
            None => None,
        }
    }
}

use Arc::*;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
/// Directed graphs.
pub struct Digraph {
    /// Number of vertices.
    size: usize,
    /// Flat matrix of arcs.
    pub edge: AntiSym<Arc>,
}

impl Digraph {
    /// Number of vertices
    pub fn size(&self) -> usize {
        self.size
    }
    
    /// Out-neigborhood of `v` in `self`.
    pub fn out_nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && self.edge.get(u, v) == BackEdge {
                res.push(u);
            }
        }
        res
    }
    /// In-neigborhood of `v` in `self`.
    pub fn in_nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && self.edge.get(u, v) == Edge {
                res.push(u);
            }
        }
        res
    }
    /// Directed graph with `n` vertices and arcs `arcs`.
    pub fn new(n: usize, arcs: &[(usize, usize)]) -> Self {
        let mut new_edge = AntiSym::new(None, n);
        for &(u, v) in arcs {
            assert!(u < n);
            assert!(v < n);
            assert!(u != v);
            assert!(new_edge.get(u, v) == None);
            new_edge.set((u, v), Edge);
        }
        Self {
            size: n,
            edge: new_edge,
        }
    }
    /// Directed graph with `n` vertices and no edge.
    pub fn empty(n: usize) -> Self {
        Self {
            size: n,
            edge: AntiSym::new(None, n),
        }
    }

    fn exist_triangle_with(&self, v: usize) -> bool {
        let out_v = self.out_nbrs(v);
        for i in self.in_nbrs(v) {
            for &o in &out_v {
                if self.edge.get(o, i) == Edge {
                    return true;
                }
            }
        }
        false
    }

    /// Directed graph obtained from `self` by adding a vertex and every edge
    /// from the rest of the graph to that vertex.
    pub fn add_sink(&self) -> Self {
        let n = self.size;
        let mut edge = self.edge.clone();
        edge.resize(n + 1, None);
        for v in 0..n {
            edge.set((v, n), Edge);
        }
        Self { edge, size: n + 1 }
    }
}

impl fmt::Display for Digraph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "(V=[{}], E={{", self.size)?;
        for u in 0..self.size {
            for v in 0..self.size {
                if v != u && self.edge.get(u, v) == Edge {
                    if self.size < 10 {
                        write!(f, " {}{}", v, u)?;
                    } else {
                        write!(f, " {}-{}", v, u)?;
                    }
                }
            }
        }
        write!(f, " }})")
    }
}

impl Canonize for Digraph {
    fn size(&self) -> usize {
        self.size
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        assert!(v < self.size);
        vec![self.out_nbrs(v), self.in_nbrs(v)]
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        self.induce(&combinatorics::invert(p))
    }
}

impl Flag for Digraph {
    fn induce(&self, p: &[usize]) -> Self {
        let k = p.len();
        let mut res = Self::empty(k);
        for u1 in 0..k {
            for u2 in 0..u1 {
                res.edge.set((u1, u2), self.edge.get(p[u1], p[u2]));
            }
        }
        res
    }

    const NAME: &'static str = "Digraph";

    fn all_flags(n: usize) -> Vec<Self> {
        if n == 0 {
            vec![Self::empty(0)]
        } else {
            unimplemented!()
        }
    }

    fn superflags(&self) -> Vec<Self> {
        let n = self.size;
        let mut res = Vec::new();
        let mut iter = Functions::new(n, 3);
        let arcs = vec![Edge, BackEdge, None];
        while let Some(f) = iter.next() {
            let mut edge = self.edge.clone();
            edge.resize(n + 1, None);
            for v in 0..n {
                edge.set((v, n), arcs[f[v]]);
            }
            res.push(Self { edge, size: n + 1 });
        }
        res
    }
}

/// Indicator for digraph without directed triangle.
#[derive(Debug, Clone, Copy)]
pub enum TriangleFree {}

impl SubFlag<Digraph> for TriangleFree {
    const SUBCLASS_NAME: &'static str = "Triangle-free digraph";
    
    fn subclass_superflags(flag: &SubClass<Digraph, Self>) -> Vec<SubClass<Digraph, Self>> {
        let n = flag.content.size;
        let mut res = Vec::new();
        let mut iter = Functions::new(n, 3);
        let arcs = vec![Edge, BackEdge, None];
        while let Some(f) = iter.next() {
            let mut edge = flag.content.edge.clone();
            edge.resize(n + 1, None);
            for v in 0..n {
                edge.set((v, n), arcs[f[v]]);
            }
            let digraph = Digraph { edge, size: n + 1 };
            if !digraph.exist_triangle_with(n) {
                res.push(digraph.into());
            }
        }
        res
    }
}

impl<A> SubClass<Digraph, A> {
    pub fn size(&self) -> usize {
        self.content.size()
    }
}
