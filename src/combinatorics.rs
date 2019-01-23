//! Basic combinatorial functions.
//!
//! In this page, `[n]` denotes the set `{0,1,...,n-1}`.
//!
//! Functions from `[n]` to `[k]` are represented by a vector of length `n`.
//! This in particular holds for permutations.
extern crate num;

use self::num::*;

/// Computes the product `start * ... * end`.
pub fn product<T>(start: T, end: T) -> T
where
    T: PrimInt,
{
    let mut res = T::one();
    for x in range(start, end + T::one()) {
        res = res * x;
    }
    res
}

/// Returns the product of integers up to the given number.
pub fn factorial<T>(n: T) -> T
where
    T: PrimInt,
{
    product(T::one(), n)
}

/// Returns the number of subsets of size `k`
/// of a set of size `n`.
pub fn binomial<T>(k: T, n: T) -> T
where
    T: PrimInt,
{
    if k < T::zero() || k > n {
        return T::zero();
    }
    product(n - k + T::one(), n) / factorial(k)
}

/// Computes the antecedents of the elements of `[n]` by `t`.
/// The output is a vector `u` of size `n`, where `u[i]` contains every `x` with `t[x]=i`.
pub fn pre_image(n: usize, t: &[usize]) -> Vec<Vec<usize>> {
    // can be improved by avoiding nested vectors
    let mut res = vec![Vec::new(); n];
    for (i, &v) in t.iter().enumerate() {
        res[v].push(i);
    }
    res
}

/// Inverts an injection from `[t.len()]` to `[n]`.
pub fn pseudo_invert(n: usize, t: &[usize]) -> Vec<Option<usize>> {
    //unsafe
    let mut res = vec![None; n];
    for (i, &v) in t.iter().enumerate() {
        debug_assert_eq!(res[v], None); // Check if t is indeed injective
        res[v] = Some(i);
    }
    res
}

/// Inverts a permutation.
pub fn invert(t: &[usize]) -> Vec<usize> {
    let n = t.len();
    let mut res = vec![n; n];
    for (i, &v) in t.iter().enumerate() {
        debug_assert_eq!(res[v], n); // Check if t is injective
        res[v] = i;
    }
    res
}

/// Assuming `t` is an injection from `[k]` to `[n]`,
/// returns the unique bijection of `[n]`
/// that is equal to `t` on `[k]` and is increasing on `{k,k+1,.., n-1}`.
pub fn permutation_of_injection(n: usize, t: &[usize]) -> Vec<usize> {
    let mut res = t.to_vec();
    let mut unassigned = vec![true; n];
    for &v in t.iter() {
        unassigned[v] = false;
    }
    for (i, &b) in unassigned.iter().enumerate() {
        if b {
            res.push(i);
        }
    }
    res
}

/*
/// Returns the array `[ p[q[0]],..., p[q[n-1]] ]`.
/// The image of `q` must be contained in `[p.len()]`.
pub fn compose(p : &[usize], q : &[usize]) -> Vec<usize> {
    debug_assert_eq!(p.len(), q.len());
    q.iter().map(|&x| p[x]).collect()
}
*/

/// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_product() {
        assert_eq!(120, product(4, 6));
        assert_eq!(17, product(17, 17));
        assert_eq!(1, product(8, 7));
    }

    #[test]
    fn unit_factorial() {
        assert_eq!(1, factorial(0));
        assert_eq!(120, factorial(5));
    }

    #[test]
    fn unit_binomial() {
        assert_eq!(35, binomial(3, 7));
        assert_eq!(0, binomial(5, 4));
        assert_eq!(0, binomial(-1, 4));
    }

    #[test]
    fn unit_pseudo_invert() {
        assert_eq!(
            pseudo_invert(5, &[3, 1, 4]),
            [None, Some(1), None, Some(0), Some(2)]
        );
    }

    #[test]
    fn unit_invert() {
        assert_eq!(invert(&[4, 0, 1, 3, 2]), [1, 2, 4, 3, 0]);
    }

    #[test]
    fn unit_permutation_of_injection() {
        assert_eq!(permutation_of_injection(6, &[2, 0, 5]), [2, 0, 5, 1, 3, 4]);
    }
}
