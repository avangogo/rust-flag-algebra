//! Undirected graphs and their implementation of `Flag`.

use crate::combinatorics;
use crate::flag::{Flag, SubFlag};
use crate::flags::common::*;
use crate::iterators;
use crate::iterators::StreamingIterator;
use canonical_form::Canonize;
use std::fmt;

/// Undirected graphs.
///
/// Implementation of graphs by upper triangular adjacency matrix,
/// reprensented as a boolean vector of length `n` choose `2`.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
pub struct Graph {
    size: usize,
    edge: SymNonRefl<bool>,
}

#[derive(Debug, Clone)]
struct EdgeIterator<'a> {
    g: &'a Graph,
    u: usize,
    v: usize,
}

impl<'a> Iterator for EdgeIterator<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.u += 1;
        if self.u >= self.v {
            self.u = 0;
            self.v += 1;
        }
        debug_assert!(self.u < self.v);
        if self.v >= self.g.size {
            None
        } else if self.g.edge(self.u, self.v) {
            Some((self.u, self.v))
        } else {
            self.next()
        }
    }
}

impl Graph {
    /// Create a graph on `n` vertices with edge set `edge`.
    /// The vertices of this graph are `0,...,n-1`.
    pub fn new(n: usize, edge: &[(usize, usize)]) -> Self {
        let mut new_edge = SymNonRefl::new(false, n);
        for &(u, v) in edge {
            assert!(u < n);
            assert!(v < n);
            new_edge[(u, v)] = true;
        }
        Self {
            size: n,
            edge: new_edge,
        }
    }
    /// Return the number of vertices in the graph
    pub fn size(&self) -> usize {
        self.size
    }

    /// Return the vector of vertices adjacent to `v`.
    pub fn nbrs(&self, v: usize) -> Vec<usize> {
        let mut res = Vec::new();
        for u in 0..self.size {
            if u != v && self.edge.get(u, v) {
                res.push(u);
            }
        }
        res
    }

    /// Create the graph on `n` vertices with no edge.  
    pub fn empty(n: usize) -> Self {
        Self {
            size: n,
            edge: SymNonRefl::new(false, n),
        }
    }

    /// Returns `true` id `uv` is an edge.
    #[inline]
    pub fn edge(&self, u: usize, v: usize) -> bool {
        u != v && self.edge[(u, v)]
    }
    /// Returns an iterator on edges.
    /// The edges are represented as couples `(u, v)` with `u < v`.
    pub fn edges(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        EdgeIterator {
            g: self,
            u: 0,
            v: 0,
        }
    }

    /// Returns `true` if the graph is connected.
    pub fn connected(&self) -> bool {
        let mut visited = vec![false; self.size];
        let mut stack = vec![0];
        while let Some(v) = stack.pop() {
            if !visited[v] {
                visited[v] = true;
                stack.append(&mut self.nbrs(v))
            }
        }
        visited.into_iter().all(|b| b)
    }
}

impl fmt::Display for Graph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "(V=[{}], E={{", self.size)?;
        for u in 0..self.size {
            for v in 0..u {
                if self.edge[(u, v)] {
                    if self.size < 10 {
                        write!(f, " {v}{u}")?
                    } else {
                        write!(f, " {v}-{u}")?
                    }
                }
            }
        }
        write!(f, " }})")
    }
}

impl Canonize for Graph {
    fn size(&self) -> usize {
        self.size
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        assert!(v < self.size);
        vec![self.nbrs(v)]
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        self.induce(&combinatorics::invert(p))
    }
}

impl Flag for Graph {
    fn induce(&self, p: &[usize]) -> Self {
        debug_assert!(p.iter().all(|&i| { i < self.size }));
        let k = p.len();
        let mut res = Self::empty(k);
        for u1 in 0..k {
            for u2 in 0..u1 {
                res.edge[(u1, u2)] = self.edge[(p[u1], p[u2])];
            }
        }
        res
    }

    const NAME: &'static str = "Graph";

    fn size_zero_flags() -> Vec<Self> {
        vec![Self::empty(0)]
    }

    fn superflags(&self) -> Vec<Self> {
        let n = self.size;
        let mut res = Vec::new();
        let mut iter = iterators::Subsets::new(n);
        while let Some(subset) = iter.next() {
            let mut edge = self.edge.clone();
            edge.resize(n + 1, false);
            for &v in subset {
                edge[(v, n)] = true;
            }
            res.push(Self { edge, size: n + 1 });
        }
        res
    }
}

// particular graphs
impl Graph {
    pub fn petersen() -> Self {
        Self::new(
            10,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 0),
                (5, 7),
                (6, 8),
                (7, 9),
                (8, 5),
                (9, 6),
                (0, 5),
                (1, 6),
                (2, 7),
                (3, 8),
                (4, 9),
            ],
        )
    }
    pub fn clique(n: usize) -> Self {
        let mut edges = Vec::new();
        for i in 0..n {
            for j in 0..i {
                edges.push((i, j))
            }
        }
        Self::new(n, &edges)
    }
    pub fn cycle(n: usize) -> Self {
        let mut edges: Vec<_> = (0..n - 1).map(|i| (i, i + 1)).collect();
        edges.push((n - 1, 0));
        Self::new(n, &edges)
    }
}

/// Connected graphs
#[derive(Debug, Clone, Copy)]
pub enum Connected {}

impl SubFlag<Graph> for Connected {
    const SUBCLASS_NAME: &'static str = "Connected graph";

    const HEREDITARY: bool = false;

    fn is_in_subclass(flag: &Graph) -> bool {
        flag.connected()
    }
}

/// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_unit() {
        let g = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
        let h = Graph::new(5, &[(3, 2), (1, 2), (3, 4), (0, 4), (1, 0)]);
        assert_eq!(g, h);
    }
    #[test]
    fn edge_iterator() {
        let g = Graph::new(5, &[(0, 1), (1, 2), (0, 4), (2, 3), (3, 4)]);
        assert_eq!(g.edges().count(), 5);
        let k3 = Graph::new(5, &[(0, 1), (1, 2), (2, 3)]);
        assert_eq!(k3.edges().count(), 3);
        let e6 = Graph::new(6, &[]);
        assert_eq!(e6.edges().count(), 0);
    }
}
