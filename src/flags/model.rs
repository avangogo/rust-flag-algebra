use crate::combinatorics::*;
use crate::common::*;
use crate::flag::Flag;
use canonical_form::Canonize;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::fmt;
use std::fmt::Display;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Serialize, Deserialize)]
pub struct Model<A> {
    size: usize,
    pub rel: A,
}

impl<A> Display for Model<A> {
    fn fmt(&self, _f: &mut fmt::Formatter) -> fmt::Result {
        unimplemented!()
    }
}

impl<R> Canonize for Model<R>
where
    R: BinRelation + Serialize + DeserializeOwned,
{
    fn size(&self) -> usize {
        self.size
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        assert!(v < self.size);
        self.rel.invariant(v)
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        self.induce(&invert(p))
    }
}

impl<R> Flag for Model<R>
where
    R: BinRelation + Serialize + DeserializeOwned,
{
    fn induce(&self, p: &[usize]) -> Self {
        Self {
            size: p.len(),
            rel: self.rel.induce(p),
        }
    }

    const NAME: &'static str = "Model";

    fn all_flags(n: usize) -> Vec<Self> {
        if n == 0 {
            vec![Self {
                size: 0,
                rel: R::empty(),
            }]
        } else {
            unimplemented!()
        }
    }

    fn superflags(&self) -> Vec<Self> {
        let size = self.size;
        self.rel
            .extensions(size)
            .into_iter()
            .map(|rel| Self {
                size: size + 1,
                rel,
            })
            .collect()
    }
}

pub type G = Model<SymNonRefl<bool>>;

/// Tests
#[cfg(test)]
mod tests {

    use super::*;

    type G = Model<SymNonRefl<bool>>;

    #[test]
    fn graph_model() {
        for (size, &nb) in [1, 1, 2, 4, 11, 34].iter().enumerate() {
            assert_eq!(G::generate(size).len(), nb);
        }
    }
}
