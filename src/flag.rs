//! Definition of the flag-able objects.

use crate::combinatorics::*;
use crate::iterators::*;
use canonical_form::*;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::cmp;
use std::collections::BTreeSet;
use std::fmt;
use std::fmt::{Debug, Display};
use std::marker::PhantomData;

/// Interface for object that can be used as flags.
pub trait Flag
where
    Self: Canonize + Debug + Display + Serialize + DeserializeOwned,
{
    /// Returns the subflag induced by the vertices in the slice `set`.
    fn induce(&self, set: &[usize]) -> Self;
    /// Returns the set of all flags of size `n`.
    ///
    /// This list can have redundancy and a priori not reduced modulo isomorphism.
    /// It is used for initialisation and does not need to be efficient.
    fn all_flags(n: usize) -> Vec<Self>;
    /// Return the list of flags of size `self.size() + 1` that contain `self`
    /// as an induced subflag.
    ///
    /// This list can have redundancy and a priori not reduced modulo isomorphism.
    fn superflags(&self) -> Vec<Self>;
    // recommended methods
    /// A name for this type of flags. For instance "Graphe".
    /// This name is used for naming the associated data subdirectory.
    fn name() -> String {
        String::from("No name")
    }

    // provided methods
    /// Return the list of flags of size `self.size() + 1` that contain `self`
    /// as an induced subflag reduced modulo isomorphism.
    fn generate_next(previous: &[Self]) -> Vec<Self> {
        let mut res: BTreeSet<Self> = BTreeSet::new();
        for g in previous {
            for h in g.superflags() {
                let _ = res.insert(canonical_form(&h));
            }
        }
        res.into_iter().collect()
    }

    /// Return the list of flags of size `n` reduced modulo isomorphism.
    fn generate(n: usize) -> Vec<Self> {
        if n == 0 {
            Self::all_flags(0)
        } else {
            Self::generate_next(&Self::generate(n - 1))
        }
    }

    /// Return the list of flags of `g_vec` that can be rooted on the
    /// flag `type_flag`.
    /// Each different way to root this flag give a different flag in the result.
    fn generate_typed_up(type_flag: &Self, g_vec: &[Self]) -> Vec<Self> {
        let n = g_vec[0].size();
        let k = type_flag.size();
        assert!(k <= n);
        let mut res: BTreeSet<Self> = BTreeSet::new();
        for g in g_vec {
            let mut iter = Injection::new(n, k);
            while let Some(selection) = iter.next() {
                if g.induce(selection) == *type_flag {
                    let p = invert(&permutation_of_injection(n, selection));
                    let _ = res.insert(canonical_form_typed(&g.apply_morphism(&p), k));
                }
            }
        }
        res.into_iter().collect()
    }

    /// Return the list of flag of size `size` rooted on `type_flag`
    /// reduced modulo (typed) isomorphism.
    fn generate_typed(type_flag: &Self, size: usize) -> Vec<Self> {
        Self::generate_typed_up(type_flag, &Self::generate(size))
    }

    /// Reorder self so that the `eta.len()` first vertices are the values
    /// of `eta` in the corresonding order.
    fn select_type(&self, eta: &[usize]) -> Self {
        let type_selector = invert(&permutation_of_injection(self.size(), eta));
        self.apply_morphism(&type_selector)
    }
}

/// A wrapper type for flags from a sub-class of flags.
pub struct SubClass<F, A> {
    /// Flag wrapped.
    pub content: F,
    phantom: PhantomData<A>,
}

impl<F, A> From<F> for SubClass<F, A> {
    #[inline]
    fn from(x: F) -> Self {
        Self {
            content: x,
            phantom: PhantomData,
        }
    }
}

impl<F: Flag, A> Clone for SubClass<F, A> {
    fn clone(&self) -> Self {
        self.content.clone().into()
    }
}
impl<F: Flag, A> Serialize for SubClass<F, A> {
    #[inline]
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.content.serialize(serializer)
    }
}
impl<'de, F: Flag, A> Deserialize<'de> for SubClass<F, A> {
    #[inline]
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        F::deserialize(deserializer).map(Into::into)
    }
}
impl<F: Flag, A> Display for SubClass<F, A> {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.content, f)
    }
}
impl<F: Flag, A> Debug for SubClass<F, A> {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Debug::fmt(&self.content, f)
    }
}
impl<F: Flag, A> PartialEq for SubClass<F, A> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.content.eq(&other.content)
    }
}
impl<F: Flag, A> Eq for SubClass<F, A> {}
impl<F: Flag, A> PartialOrd for SubClass<F, A> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<F: Flag, A> Ord for SubClass<F, A> {
    #[inline]
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.content.cmp(&other.content)
    }
}

/// Mechanism for defining a subclass of a flag class.
pub trait SubFlag<F>
where
    F: Flag,
    Self: Sized,
{
    fn subclass_superflags(flag: &SubClass<F, Self>) -> Vec<SubClass<F, Self>>;
    fn subclass_name() -> String;
}

impl<F, A> Canonize for SubClass<F, A>
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
    fn apply_morphism(&self, p: &[usize]) -> Self {
        self.content.apply_morphism(p).into()
    }
}

impl<F, A> Flag for SubClass<F, A>
where
    A: SubFlag<F>,
    F: Flag,
{
    fn name() -> String {
        A::subclass_name()
    }
    fn superflags(&self) -> Vec<Self> {
        A::subclass_superflags(self)
    }
    // inherited methods
    fn induce(&self, p: &[usize]) -> Self {
        self.content.induce(p).into()
    }
    fn all_flags(n: usize) -> Vec<Self> {
        F::all_flags(n).into_iter().map(Into::into).collect()
    }
}
