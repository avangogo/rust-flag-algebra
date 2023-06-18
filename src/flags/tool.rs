use crate::flag::Flag;
use canonical_form::Canonize;
use std::cmp::*;
use std::fmt;
use std::fmt::{Debug, Display};

/// A type for colored flags
#[derive(Clone, PartialOrd, PartialEq, Ord, Eq, Debug, Hash, Serialize, Deserialize)]
pub struct Colored<F, const N: u8> {
    pub content: F,
    pub color: Vec<u8>,
}

impl<F, const N: u8> Colored<F, N>
where
    F: Canonize,
{
    pub fn new(content: F, color: Vec<u8>) -> Self {
        assert_eq!(content.size(), color.len());
        for &c in &color {
            assert!(c < N)
        }
        Self { content, color }
    }
}

impl<F: Flag, const N: u8> Display for Colored<F, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.content, f)?;
        Debug::fmt(&self.color, f)
    }
}

impl<F, const N: u8> Flag for Colored<F, N>
where
    F: Flag,
{
    /// Returns the subflag induced by the vertices in the slice `set`.
    fn induce(&self, set: &[usize]) -> Self {
        Self {
            content: self.content.induce(set),
            color: set.iter().map(|u| self.color[*u]).collect(),
        }
    }
    fn size_zero_flags() -> Vec<Self> {
        F::size_zero_flags()
            .into_iter()
            .map(|flag: F| Self {
                content: flag,
                color: Vec::new(),
            })
            .collect()
    }
    fn superflags(&self) -> Vec<Self> {
        let mut res = Vec::new();
        for flag in self.content.superflags().into_iter() {
            for c in 0..N {
                let mut color = self.color.clone();
                color.push(c);
                res.push(Self {
                    content: flag.clone(),
                    color,
                })
            }
        }
        res
    }

    // FIXME: Incorrect name because of limitation of consts in Rusts
    // Should be 'fmt("{}-colored {}", N, F::NAME'
    const NAME: &'static str = "FIXME";
    const HEREDITARY: bool = F::HEREDITARY;
}

impl<F, const N: u8> Canonize for Colored<F, N>
where
    F: Flag,
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
                    col[i] = col[i] * (N as u64) + (*c as u64);
                }
                Some(col)
            }
        }
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        let mut color: Vec<u8> = vec![N; p.len()];
        for (i, &pi) in p.iter().enumerate() {
            color[pi] = self.color[i]
        }
        Self {
            content: self.content.apply_morphism(p),
            color,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::Graph;
    use canonical_form::Canonize;
    #[test]
    fn test_colored() {
        type G2 = Colored<Graph, 2>;
        let p3 = Graph::new(3, &[(0, 1), (1, 2)]);
        let g: G2 = Colored {
            content: p3.clone(),
            color: vec![0, 0, 1],
        };
        let h: G2 = Colored {
            content: p3,
            color: vec![1, 0, 0],
        };
        assert_eq!(g.canonical(), h.canonical());
        type G1 = Colored<Graph, 1>;
        assert_eq!(G1::generate(5).len(), 34);
        type G5 = Colored<Graph, 5>;
        assert_eq!(G5::generate(2).len(), 30);
        assert_eq!(G2::generate(3).len(), 20);
    }
}
