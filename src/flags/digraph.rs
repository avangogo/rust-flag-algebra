//! `[OrientedGraph]`s and their flag implementation.
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
    /// Arcs in both directions
    Reciprocal,
}

impl Neg for Arc {
    type Output = Self;

    fn neg(self) -> Self {
        match self {
            Edge => BackEdge,
            BackEdge => Edge,
            e => e,
        }
    }
}

use Arc::*;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
/// Loopless directed graphs.
///
/// The possible relations between two vertices is described by the type [`Arc`].
pub struct DirectedGraph {
    /// Number of vertices.
    size: usize,
    /// Flat matrix of arcs.
    edge: AntiSym<Arc>,
}

impl DirectedGraph {
    /// Create a directed graph with `n` vertices and arcs in `arcs`.
    ///
    /// # Panics
    /// * If `arcs` contains vertices not in `{0, ..., n-1}`
    /// * If a loop `(u, u)` is provided
    /// * If some arc is provided twice
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    ///
    /// // Oriented path with 3 vertices `0 -> 1 -> 2`
    /// let p3 = OrientedGraph::new(3, [(0, 1), (1, 2)]);
    /// ```
    pub fn new<I>(n: usize, arcs: I) -> Self
    where
        I: IntoIterator<Item = (usize, usize)>,
    {
        let mut new_edge = AntiSym::new(None, n);
        for (u, v) in arcs {
            check_arc((u, v), n);
            let new_arc = match new_edge.get(u, v) {
                None => Edge,
                BackEdge => Reciprocal,
                _ => panic!("Arc ({u}, {v}) specified twice"),
            };
            new_edge.set((u, v), new_arc);
        }
        Self {
            size: n,
            edge: new_edge,
        }
    }
    /// Number of vertices
    pub fn size(&self) -> usize {
        self.size
    }
    /// Out-neigborhood of `v`.
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    /// let p3 = OrientedGraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.out_nbrs(1), vec![2]);
    /// ```
    pub fn out_nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && matches!(self.edge.get(u, v), BackEdge | Reciprocal) {
                res.push(u);
            }
        }
        res
    }
    /// In-neigborhood of `v`.
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    /// let p3 = OrientedGraph::new(3, [(0, 1), (1, 2)]);
    /// assert_eq!(p3.in_nbrs(1), vec![0]);
    /// ```
    pub fn in_nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && matches!(self.edge.get(u, v), Edge | Reciprocal) {
                res.push(u);
            }
        }
        res
    }
    /// Oriented relation between `u` to `v`.
    /// ```
    /// use flag_algebra::flags::{DirectedGraph, Arc};
    /// let g = DirectedGraph::new(3, [(0, 1), (1, 2), (2, 1)]);
    /// assert_eq!(g.arc(0, 1), Arc::Edge);
    /// assert_eq!(g.arc(1, 0), Arc::BackEdge);
    /// assert_eq!(g.arc(0, 2), Arc::None);
    /// assert_eq!(g.arc(1, 2), Arc::Reciprocal);
    /// ```
    pub fn arc(&self, u: usize, v: usize) -> Arc {
        check_arc((u, v), self.size);
        self.edge.get(u, v)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
/// Loopless oriented graphs.
///
/// Every pair of vertices has at most one edge.
pub struct OrientedGraph(DirectedGraph);

impl OrientedGraph {
    pub fn as_directed(&self) -> &DirectedGraph {
        &self.0
    }
    /// Create an oriented graph with `n` vertices and arcs in `arcs`.
    ///
    /// # Panics
    /// * If `arcs` contains vertices not in `{0, ..., n-1}`
    /// * If a loop `(u, u)` is provided
    /// * If some arc is provided twice
    /// * If both arcs `(u, v)` and `(v, u)` are provided
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    ///
    /// // Oriented path with 3 vertices `0 -> 1 -> 2`
    /// let p3 = OrientedGraph::new(3, [(0, 1), (1, 2)]);
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
        Self(DirectedGraph {
            size: n,
            edge: new_edge,
        })
    }
    /// Oriented graph with `n` vertices and no edge.
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    /// assert_eq!(OrientedGraph::empty(4), OrientedGraph::new(4, []));
    /// ```
    pub fn empty(n: usize) -> Self {
        Self::new(n, [])
    }
    /// Number of vertices
    pub fn size(&self) -> usize {
        self.0.size()
    }
    /// Out-neigborhood of `v`.
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    /// let p3 = OrientedGraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.out_nbrs(1), vec![2]);
    /// ```
    pub fn out_nbrs(&self, v: usize) -> Vec<usize> {
        self.0.out_nbrs(v)
    }
    /// In-neigborhood of `v`.
    /// ```
    /// use flag_algebra::flags::OrientedGraph;
    /// let p3 = OrientedGraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.in_nbrs(1), vec![0]);
    /// ```
    pub fn in_nbrs(&self, v: usize) -> Vec<usize> {
        self.0.in_nbrs(v)
    }
    /// Oriented relation between `u` to `v`.
    /// ```
    /// use flag_algebra::flags::{OrientedGraph, Arc};
    /// let p3 = OrientedGraph::new(3, [(0,1), (1,2)]);
    /// assert_eq!(p3.arc(0, 1), Arc::Edge);
    /// assert_eq!(p3.arc(1, 0), Arc::BackEdge);
    /// assert_eq!(p3.arc(0, 2), Arc::None);
    /// ```
    pub fn arc(&self, u: usize, v: usize) -> Arc {
        self.0.arc(u, v)
    }
    fn is_triangle_free(&self) -> bool {
        for u in 0..self.size() {
            // Assume u is the largest vertex
            for v in 0..u {
                if self.arc(v, u) == Edge {
                    for w in 0..u {
                        if self.arc(u, w) == Edge && self.arc(w, v) == Edge {
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
        let n = self.size();
        let mut edge = self.0.edge.clone();
        edge.resize(n + 1, None);
        for v in 0..n {
            edge.set((v, n), Edge);
        }
        Self(DirectedGraph { edge, size: n + 1 })
    }
}

impl fmt::Display for DirectedGraph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "(V=[{}], E={{", self.size)?;
        for u in 0..self.size {
            for v in 0..self.size {
                if v != u {
                    match self.edge.get(u, v) {
                        Edge => write!(f, " {u}->{v}")?,
                        Reciprocal if u < v => write!(f, " {u}<->{v}")?,
                        _ => (),
                    }
                }
            }
        }
        write!(f, " }})")
    }
}

impl fmt::Display for OrientedGraph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.fmt(f)
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
    assert!(
        u != v,
        "Invalid arc ({u}, {v}): OrientedGraphs have no loop"
    );
}

impl Canonize for DirectedGraph {
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

impl Canonize for OrientedGraph {
    fn size(&self) -> usize {
        self.size()
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        self.0.invariant_neighborhood(v)
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        Self(self.0.apply_morphism(p))
    }
}

impl Flag for DirectedGraph {
    fn induce(&self, p: &[usize]) -> Self {
        let k = p.len();
        let mut res = Self::new(k, []);
        for u1 in 0..k {
            for u2 in 0..u1 {
                res.edge.set((u1, u2), self.edge.get(p[u1], p[u2]));
            }
        }
        res
    }

    const NAME: &'static str = "DirectedGraph";

    fn size_zero_flags() -> Vec<Self> {
        vec![Self::new(0, [])]
    }

    fn superflags(&self) -> Vec<Self> {
        let mut res = Vec::new();
        let mut iter = Functions::new(self.size, 4);
        let arcs = [Edge, BackEdge, None, Reciprocal];
        while let Some(f) = iter.next() {
            res.push(extend(self, |v| arcs[f[v]]));
        }
        res
    }
}

fn extend<F: Fn(usize) -> Arc>(g: &DirectedGraph, f: F) -> DirectedGraph {
    let mut edge = g.edge.clone();
    let n = g.size;
    edge.resize(n + 1, None);
    for v in 0..n {
        edge.set((v, n), f(v));
    }
    DirectedGraph { size: n + 1, edge }
}

impl Flag for OrientedGraph {
    fn induce(&self, p: &[usize]) -> Self {
        Self(self.0.induce(p))
    }

    const NAME: &'static str = "OrientedGraph";

    fn size_zero_flags() -> Vec<Self> {
        vec![Self::empty(0)]
    }

    fn superflags(&self) -> Vec<Self> {
        let mut res = Vec::new();
        let mut iter = Functions::new(self.size(), 3);
        let arcs = [Edge, BackEdge, None];
        while let Some(f) = iter.next() {
            res.push(Self(extend(&self.0, |v| arcs[f[v]])));
        }
        res
    }
}

/// Indicator for oriented graph without directed triangle.
///
/// Makes `SubClass<OrientedGraph, TriangleFree>` usable as flags.
#[derive(Debug, Clone, Copy)]
pub enum TriangleFree {}

impl SubFlag<OrientedGraph> for TriangleFree {
    const SUBCLASS_NAME: &'static str = "Triangle-free oriented graph";

    fn is_in_subclass(flag: &OrientedGraph) -> bool {
        flag.is_triangle_free()
    }
}

impl<A> SubClass<OrientedGraph, A> {
    pub fn size(&self) -> usize {
        self.content.size()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display() {
        let g = DirectedGraph::new(3, [(0, 1), (1, 2), (2, 1)]);
        assert_eq!(format!("{g}"), "(V=[3], E={ 0->1 1<->2 })");
    }

    #[test]
    fn test_display_oriented() {
        let g = OrientedGraph::new(2, [(0, 1)]);
        assert_eq!(format!("{g}"), "(V=[2], E={ 0->1 })");
    }
}
