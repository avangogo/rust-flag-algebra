//! `[Digraph]`s and their flag implementation.
use crate::combinatorics;
use crate::flag::{Flag, SubClass, SubFlag};
use crate::flags::common::*;
use crate::iterators::{Functions, StreamingIterator};
use canonical_form::Canonize;
use std::fmt;
use std::ops::Neg;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy, Serialize, Deserialize)]
/// The arc between two nodes of a directed graphs.
#[derive(Default)]
pub enum Arc {
    /// No arc.
    #[default]
    None,
    /// Arc from the first to the second vertex considered.
    Edge,
    /// Arc from the second to the first vertex considered.
    BackEdge,
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
/// Loopless oriented graphs.
///
/// Every pair of vertices has at most one edge.
/// The possible relations between two vertices is described by the type [`Arc`].
pub struct Digraph {
    /// Number of vertices.
    size: usize,
    /// Flat matrix of arcs.
    edge: AntiSym<Arc>,
}

impl Digraph {
    /// Create an oriented graph with `n` vertices and arcs in `arcs`.
    ///
    /// # Panics
    /// * If `arcs` contains vertices not in `{0, ..., n-1}`
    /// * If a loop `(u, u)` is provided
    /// * If some arc is provided twice
    /// * If both arcs `(u, v)` and `(v, u)` are provided
    /// ```
    /// use flag_algebra::flags::Digraph;
    ///
    /// // Oriented path with 3 vertices `0 -> 1 -> 2`
    /// let p3 = Digraph::new(3, [(0, 1), (1, 2)]);
    /// ```
    pub fn new<I>(n: usize, arcs: I) -> Self
    where
        I: IntoIterator<Item = (usize, usize)>,
    {
        let mut new_edge = AntiSym::new(None, n);
        for (u, v) in arcs {
            check_arc((u, v), n);
            assert!(
                new_edge.get(u, v) == None,
                "Pair {{{u}, {v}}} specified twice"
            );
            new_edge.set((u, v), Edge);
        }
        Self {
            size: n,
            edge: new_edge,
        }
    }
    /// Oriented graph with `n` vertices and no edge.
    /// ```
    /// use flag_algebra::flags::Digraph;
    /// assert_eq!(Digraph::empty(4), Digraph::new(4, []));
    /// ```
    pub fn empty(n: usize) -> Self {
        Self {
            size: n,
            edge: AntiSym::new(None, n),
        }
    }
    /// Number of vertices
    pub fn size(&self) -> usize {
        self.size
    }
    /// Out-neigborhood of `v`.
    /// ```
    /// use flag_algebra::flags::Digraph;
    /// let p3 = Digraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.out_nbrs(1), vec![2]);
    /// ```
    pub fn out_nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && self.edge.get(u, v) == BackEdge {
                res.push(u);
            }
        }
        res
    }
    /// In-neigborhood of `v`.
    /// ```
    /// use flag_algebra::flags::Digraph;
    /// let p3 = Digraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.in_nbrs(1), vec![0]);
    /// ```
    pub fn in_nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && self.edge.get(u, v) == Edge {
                res.push(u);
            }
        }
        res
    }
    /// Oriented relation between `u` to `v`.
    /// ```
    /// use flag_algebra::flags::{Digraph, Arc};
    /// let p3 = Digraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.arc(0, 1), Arc::Edge);
    /// assert_eq!(p3.arc(1, 0), Arc::BackEdge);
    /// assert_eq!(p3.arc(0, 2), Arc::None);
    /// ```
    pub fn arc(&self, u: usize, v: usize) -> Arc {
        check_arc((u, v), self.size);
        self.edge.get(u, v)
    }
    fn is_triangle_free(&self) -> bool {
        for u in 0..self.size {
            // Assume u is the largest vertex
            for v in 0..u {
                if self.edge.get(v, u) == Edge {
                    for w in 0..u {
                        if self.edge.get(u, w) == Edge && self.edge.get(w, v) == Edge {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Oriented graph obtained from `self` by adding a vertex and every edge
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

fn check_vertex(u: usize, graph_size: usize) {
    assert!(
        u < graph_size,
        "Invalid vertex {u}: the vertex set is {{0, ..., {}}}",
        graph_size - 1
    );
}

fn check_arc((u, v): (usize, usize), graph_size: usize) {
    check_vertex(u, graph_size);
    check_vertex(v, graph_size);
    assert!(u != v, "Invalid arc ({u}, {v}): Digraph have no loop");
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

    fn size_zero_flags() -> Vec<Self> {
        vec![Self::empty(0)]
    }

    fn superflags(&self) -> Vec<Self> {
        let n = self.size;
        let mut res = Vec::new();
        let mut iter = Functions::new(n, 3);
        let arcs = [Edge, BackEdge, None];
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
///
/// Makes `SubClass<Digraph, TriangleFree>` usable as flags.
#[derive(Debug, Clone, Copy)]
pub enum TriangleFree {}

impl SubFlag<Digraph> for TriangleFree {
    const SUBCLASS_NAME: &'static str = "Triangle-free digraph";

    fn is_in_subclass(flag: &Digraph) -> bool {
        flag.is_triangle_free()
    }
}

impl<A> SubClass<Digraph, A> {
    pub fn size(&self) -> usize {
        self.content.size()
    }
}
