//use crate::common::*;
use crate::flag::Flag; //, SubClass, SubFlag};
use canonical_form::Canonize;
use core::marker::PhantomData;
use serde::de::DeserializeOwned;
use serde::Serialize;
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

impl<F: Flag, N> Display for Colored<F, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.content, f)?;
        Debug::fmt(&self.color, f)
    }
}

impl<F, N> Flag for Colored<F, N>
where
    F: Flag,
    N: Unsigned + Clone + Eq + Ord + Debug + Serialize + DeserializeOwned,
{
    /// Returns the subflag induced by the vertices in the slice `set`.
    fn induce(&self, set: &[usize]) -> Self {
        Self {
            content: self.content.induce(set),
            color: set.into_iter().map(|u| self.color[*u]).collect(),
            phantom: PhantomData,
        }
    }
    fn all_flags(n: usize) -> Vec<Self> {
        if n == 0 {
            F::all_flags(0)
                .into_iter()
                .map(|flag: F| Self {
                    content: flag,
                    color: Vec::new(),
                    phantom: PhantomData,
                })
                .collect()
        } else {
            unimplemented!()
        }
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
    N: Unsigned + Clone + Eq + Ord + Debug + Serialize + DeserializeOwned,
{
    #[inline]
    fn size(&self) -> usize {
        self.content.size()
    }
    fn invariant_neighborhood(&self, v: usize) -> Vec<Vec<usize>> {
        self.content.invariant_neighborhood(v)
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        Self {
            content: self.content.apply_morphism(p),
            color: p.iter().map(|i| self.color[*i]).collect(),
            phantom: PhantomData,
        }
    }
}
