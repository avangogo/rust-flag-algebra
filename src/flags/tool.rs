use crate::flag::Flag;
use canonical_form::Canonize;
use core::marker::PhantomData;
use std::cmp::*;
use std::fmt;
use std::fmt::{Debug, Display};
use typenum::Unsigned;

/// A type for colored flags
#[derive(Clone, PartialOrd, PartialEq, Ord, Eq, Debug, Hash, Serialize, Deserialize)]
pub struct Colored<F, N> {
    pub content: F,
    pub color: Vec<u8>,
    phantom: PhantomData<N>,
}

impl<F, N> Colored<F, N>
where
    N: Unsigned,
    F: Canonize,
{
    pub fn new(content: F, color: Vec<u8>) -> Self {
        assert_eq!(content.size(), color.len());
        for &c in &color {
            assert!(c < N::U8)
        }
        Self {
            content,
            color,
            phantom: PhantomData,
        }
    }
}

impl<F: Flag, N> Display for Colored<F, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.content, f)?;
        Debug::fmt(&self.color, f)
    }
}

impl<F, N> Flag for Colored<F, N>
where
    F: Flag,
    N: Unsigned + Clone + Eq + Ord + Debug,
{
    /// Returns the subflag induced by the vertices in the slice `set`.
    fn induce(&self, set: &[usize]) -> Self {
        Self {
            content: self.content.induce(set),
            color: set.iter().map(|u| self.color[*u]).collect(),
            phantom: PhantomData,
        }
    }
    fn size_zero_flags() -> Vec<Self> {
        F::size_zero_flags()
            .into_iter()
            .map(|flag: F| Self {
                content: flag,
                color: Vec::new(),
                phantom: PhantomData,
            })
            .collect()
    }
    fn superflags(&self) -> Vec<Self> {
        let mut res = Vec::new();
        for flag in self.content.superflags().into_iter() {
            for c in 0..N::U8 {
                let mut color = self.color.clone();
                color.push(c);
                res.push(Self {
                    content: flag.clone(),
                    color,
                    phantom: PhantomData,
                })
            }
        }
        res
    }

    const NAME: &'static str = "TODO";
    const HEREDITARY: bool = F::HEREDITARY;
}

impl<F, N> Canonize for Colored<F, N>
where
    F: Flag,
    N: Unsigned + Clone + Eq + Ord,
{
    #[inline]
    fn size(&self) -> usize {
        self.content.size()
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        self.content.invariant_neighborhood(v)
    }
    fn invariant_coloring(&self) -> Option<Vec<u64>> {
        match self.content.invariant_coloring() {
            None => Some(self.color.iter().map(|&c| c as u64).collect()),
            Some(mut col) => {
                for (i, c) in self.color.iter().enumerate() {
                    col[i] = col[i] * N::U64 + (*c as u64);
                }
                Some(col)
            }
        }
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        let mut color: Vec<u8> = vec![N::U8; p.len()];
        for (i, &pi) in p.iter().enumerate() {
            color[pi] = self.color[i]
        }
        Self {
            content: self.content.apply_morphism(p),
            color,
            phantom: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::Graph;
    use canonical_form::Canonize;
    use typenum::*;
    #[test]
    fn test_colored() {
        type G2 = Colored<Graph, U2>;
        let p3 = Graph::new(3, &[(0, 1), (1, 2)]);
        let g: G2 = Colored {
            content: p3.clone(),
            color: vec![0, 0, 1],
            phantom: PhantomData,
        };
        let h: G2 = Colored {
            content: p3,
            color: vec![1, 0, 0],
            phantom: PhantomData,
        };
        assert_eq!(g.canonical(), h.canonical());
        type G1 = Colored<Graph, U1>;
        assert_eq!(G1::generate(5).len(), 34);
        type G5 = Colored<Graph, U5>;
        assert_eq!(G5::generate(2).len(), 30);
        assert_eq!(G2::generate(3).len(), 20);
    }
}
