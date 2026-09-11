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
        for flag in self.content.superflags() {
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

    fn name() -> String {
        format!("{N}-colored {}", F::name())
    }
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
    fn invariant_neighborhood(&self, v: usize) -> impl Iterator<Item = (usize, u64)> {
        self.content.invariant_neighborhood(v)
    }
    fn invariant_color(&self, v: usize) -> u64 {
        self.color[v] as u64 + (N as u64).wrapping_mul(self.content.invariant_color(v))
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
    use crate::operator::{Basis, Savable};
    use canonical_form::Canonize;

    /// Every `Colored` class used to report the name "FIXME", so all of them
    /// shared one memoization directory and silently served each other's flag
    /// lists. Names must distinguish both the colour count and the base class.
    #[test]
    fn colored_names_are_distinct() {
        assert_eq!(<Colored<Graph, 2> as Flag>::name(), "2-colored Graph");
        assert_eq!(<Colored<Graph, 3> as Flag>::name(), "3-colored Graph");
        assert_eq!(
            <Colored<Colored<Graph, 2>, 3> as Flag>::name(),
            "3-colored 2-colored Graph"
        );
    }

    /// The cached bases must agree with `generate`, which never reads the cache.
    #[test]
    fn colored_bases_are_not_shared() {
        for size in 0..3 {
            assert_eq!(
                Basis::<Colored<Graph, 2>>::new(size).get(),
                <Colored<Graph, 2> as Flag>::generate(size)
            );
            assert_eq!(
                Basis::<Colored<Graph, 3>>::new(size).get(),
                <Colored<Graph, 3> as Flag>::generate(size)
            );
        }
    }
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
