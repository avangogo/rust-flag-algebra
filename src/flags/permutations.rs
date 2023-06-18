//! Permutations.

use crate::combinatorics;
use crate::common::*;
use crate::flag::{Flag, SubClass};
use crate::flags::Colored;
use crate::iterators;
use crate::iterators::StreamingIterator;
use canonical_form::Canonize;
use std::fmt;
use std::fmt::Debug;

/// Permutations.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
pub struct Permutation(Vec<usize>);

fn is_permutation(p: &[u8]) -> bool {
    let n = p.len();
    let seen = vec![false; n];
    for i_u8 in p {
        let i = i_u8 as usize; 
        if i >= n || seen[i] {
            return false
        }
        seen[i] = true;
    }
    true
}

/*
impl From<&str> for Permutation {
    fn from(s: &str) -> Self {
        let data: Vec<u8> = s.iter().map()
    }
}
*/

impl Permutation {
    /// Create a new permutation.
    pub fn new(data: &[u8]) -> Self {
        if !is_permutation(data) {
            panic!("Not a permutation: {:?}", data)
        }
        Self(data.to_vec())
    }    
}

impl<E> fmt::Display for Permutations {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in self.data {
            write!(f, "{}", i)
        }
    } 
}

// Canonization is trivial for permutations
impl Canonize for Permutation
{
    fn size(&self) -> usize {
        self.0.len()
    }
    fn invariant_neighborhood(&self, _: usize) -> Vec<Vec<usize>> {
        Vec::new()
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

    // Should be something like `concat!("CGraph_", E::USIZE)` which currently does not work
    // Waiting for further support of `const fn`
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
    }
}
