//! Definition of the flag-able objects.

use crate::combinatorics::*;
use crate::iterators::*;
use canonical_form::Canonize;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::cmp;
use std::collections::BTreeSet;
use std::fmt;
use std::fmt::{Debug, Display};
use std::marker::PhantomData;

/// Trait for combinatorial objects that can be used as flags.
///
/// Flags must implement the `Canonize` trait
/// from the `canonical_form` crate, that allows
/// reduction modulo isomorphism.
pub trait Flag
where
    Self: Canonize + Debug + Display + Serialize + DeserializeOwned,
{
    /// Returns the subflag induced by the vertices in the slice `set`.
    fn induce(&self, set: &[usize]) -> Self;
    /// Returns the set of all flags of size 0.
    ///
    /// This list is used for initialisation in the flag construction.
    fn size_zero_flags() -> Vec<Self>;
    /// Return the list of flags of size `self.size() + 1` that contain `self`
    /// as an induced subflag.
    ///
    /// This list can have redundancy and is a priori not reduced modulo isomorphism.    
    fn superflags(&self) -> Vec<Self>;
    /// A unique name for this type of flags. For instance "Graph".
    /// This nameis used for naming the associated data subdirectory.
    const NAME: &'static str;

    // caracteristic
    /// Setting this parameter to `false` deactivate checks that induced subflags exists.
    /// Must be `true` in every classic case.
    const HEREDITARY: bool = true;

    // provided methods
    /// Return the list of flags of size `self.size() + 1` that contain `self`
    /// as an induced subflag reduced modulo isomorphism.
    fn generate_next(previous: &[Self]) -> Vec<Self> {
        let mut res: BTreeSet<Self> = BTreeSet::new();
        for g in previous {
            for h in g.superflags() {
                let _ = res.insert(h.canonical());
            }
        }
        res.into_iter().collect()
    }

    /// Return the list of flags of size `n` reduced modulo isomorphism.
    fn generate(n: usize) -> Vec<Self> {
        if n == 0 {
            Self::size_zero_flags()
        } else {
            Self::generate_next(&Self::generate(n - 1))
        }
    }

    /// Return the list of flags of `g_vec` that can be rooted on the
    /// flag `type_flag`.
    /// Each different way to root this flag give a different flag in the result.
    fn generate_typed_up(type_flag: &Self, g_vec: &[Self]) -> Vec<Self> {
        assert!(!g_vec.is_empty());
        let n = g_vec[0].size();
        let k = type_flag.size();
        assert!(k <= n);
        let mut res: BTreeSet<Self> = BTreeSet::new();
        for g in g_vec {
            // For every subset of size k
            let mut iter = Choose::new(n, k);
            while let Some(pre_selection) = iter.next() {
                // If this subset induces the right type, then test all permutations
                if &g.induce(pre_selection).canonical() == type_flag {
                    let mut iter2 = Injection::permutation(k);
                    while let Some(select2) = iter2.next() {
                        let selection = &compose(pre_selection, select2);
                        let h = g.induce(selection);
                        if &h == type_flag {
                            let p = invert(&permutation_of_injection(n, selection));
                            let _ = res.insert(g.apply_morphism(&p).canonical_typed(k));
                        }
                    }
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
    /// of `eta` in the corresponding order.
    fn select_type(&self, eta: &[usize]) -> Self {
        let type_selector = invert(&permutation_of_injection(self.size(), eta));
        self.apply_morphism(&type_selector)
    }
}

/// A wrapper type for flags from a sub-class of flags.
///
/// This structure is meant to be used with the [`SubFlag`](trait.SubFlag.html) trait.
/// The second type parameter serves as an identifier for the subclass
/// and should implement `SubFlag<F>`.
///
/// See the [`SubFlag`](trait.SubFlag.html) page for an example.
pub struct SubClass<F, A> {
    /// Type of flag wrapped.
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
    #[inline]
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

/// Trait for defining a class of flag as a subset of an already defined class of flag.
///
/// The dimension of the related flag algebra will be equal to the size of the subclass.
///
/// Flags of a sublass are generated by filtering flags with the
/// `is_in_subclass` method when extending
/// Flags of size `n` from flags of size `n-1` (in the subclass).
/// Since  filtering happens at each step, the exhaustive list of flags of
/// the general class is never generated.
///
/// # Example
/// ```
/// // This shows how to define triangle-free graphs based on graphs.
/// use flag_algebra::*;
/// use flag_algebra::flags::Graph;
///
/// // We first define a (zero-sized) type identifying the subclass
/// enum TriangleFree {}
///
/// // We want `SubClass<Graph, TriangleFree>` to be a subclass of `Graph`.
/// // This is done by implementing `SubFlag<Graph>` for `TriangleFree`.
/// impl SubFlag<Graph> for TriangleFree {
///     const SUBCLASS_NAME: &'static str = "Triangle-free graph for the example";
///
///     // Compute if the graph is triangle-free.
///     fn is_in_subclass(g: &Graph) -> bool {
///         for (u, v) in g.edges() {
///             for w in 0..u {
///                 if g.edge(u, w) && g.edge(v, w) {
///                     return false // Found a triangle
///                 }
///             }
///         }
///         true
///     }
/// }
///
/// // We can now use `SubClass<Graph, TriangleFree>` as flags for triangle-free graphs.
/// type F = SubClass<Graph, TriangleFree>;
/// let basis: Basis<F> = Basis::new(4);
/// assert_eq!(basis.get().len(), 7); // There are 7 triangle-free graphs of size 4
/// ```
pub trait SubFlag<F>
where
    F: Flag,
    Self: Sized,
{
    /// Identifier function for the subclass.
    fn is_in_subclass(flag: &F) -> bool;

    /// Unique name for the subclass.
    /// This is used for naming the memoization directory.
    const SUBCLASS_NAME: &'static str;

    const HEREDITARY: bool = F::HEREDITARY;
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
    const NAME: &'static str = A::SUBCLASS_NAME;

    const HEREDITARY: bool = A::HEREDITARY;

    fn superflags(&self) -> Vec<Self> {
        let mut res = Vec::new();
        for flag in self.content.superflags() {
            if A::is_in_subclass(&flag) {
                res.push(flag.into())
            }
        }
        res
    }
    // inherited methods
    fn induce(&self, p: &[usize]) -> Self {
        self.content.induce(p).into()
    }
    fn size_zero_flags() -> Vec<Self> {
        F::size_zero_flags().into_iter().map(Into::into).collect()
    }
}
