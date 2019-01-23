//! Computing coefficients of a flag algebra operator.

use crate::flag::Flag;
use crate::iterators::*;
use canonical_form::*;
use sprs::{CsMat, TriMatI, CSC, CSR};

fn induce_and_reduce<F: Flag>(type_size: usize, f: &F, subset: &[usize]) -> F {
    canonical_form_typed(&f.induce(subset), type_size)
}

/// Returns the number of induced subflags of `g` isomorphic to `f`,
/// considered with type of size `type_size`.
pub fn count_subflags<F: Flag>(type_size: usize, f: &F, g: &F) -> u64 {
    // f must be in normal form
    assert_eq!(*f, canonical_form_typed(f, type_size));
    let k = f.size();
    let n = g.size();
    assert!(type_size <= k && k <= n);
    let mut count: u64 = 0;
    let mut iter = Choose::with_fixed_part(n, k, type_size);
    while let Some(subset) = iter.next() {
        if induce_and_reduce(type_size, g, subset) == *f {
            count += 1;
        }
    }
    count
}

/// Returns the number of split partitions of `g`
/// with intersection `[type_size]`
/// inducing `f1` and `f2`.
pub fn count_split<F: Flag>(type_size: usize, f1: &F, f2: &F, g: &F) -> u64 {
    assert_eq!(*f1, canonical_form_typed(f1, type_size));
    assert_eq!(*f2, canonical_form_typed(f2, type_size));
    let k1 = f1.size();
    let k2 = f2.size();
    let n = g.size();
    assert!(type_size <= k1 && k1 <= n);
    assert_eq!(n, k1 + k2 - type_size);
    let mut count: u64 = 0;
    let mut iter = Split::with_fixed_part(n, k1, type_size);
    while let Some((subset1, subset2)) = iter.next() {
        if *f1 == induce_and_reduce(type_size, g, subset1)
            && *f2 == induce_and_reduce(type_size, g, subset2)
        {
            count += 1;
        }
    }
    count
}

/// Returns a matrix `m`
/// where `m[i,j]` is equal to
/// `count_subflags(type_size, f_vec[i], g_vec[j])`.
///
/// This more efficient than iterating `count_subflags`.
///
/// The output is a sparse matrix in CSR form.
pub fn count_subflag_tabulate<F: Flag>(type_size: usize, f_vec: &[F], g_vec: &[F]) -> CsMat<u64> {
    let k = f_vec[0].size();
    let n = g_vec[0].size();
    assert!(type_size <= k && k <= n);
    //
    let mut res = CsMat::empty(CSR, f_vec.len());
    for g in g_vec {
        let mut column = vec![0; f_vec.len()];
        //
        let mut iter = Choose::with_fixed_part(n, k, type_size);
        while let Some(subset) = iter.next() {
            let f0 = induce_and_reduce(type_size, g, subset);
            let i = f_vec
                .binary_search(&f0)
                .unwrap_or_else(|_| panic!("Flag not found: {:?}", &f0));
            column[i] += 1;
        }
        res = res.append_outer(&column)
    }
    res
}

/// Returns a vector `m` of length `g_vec.len()` of
/// `f1_vec.len()`x`f2_vec.len()` matrices
/// such that
/// `m[i][j,k]` is equal to
/// `count_split(type_size, f1_vec[j], f2_vec[k], g_vec[i])`.
///
/// This more efficient than iterating `count_split`.
///
/// The matrices are in CSR form.
// So that Trans(f1) * M * f2 is correct
// reminder: CSR: outer = rows, inner = cols
// so it is a outer x inner matrix (rows x cols)
pub fn count_split_tabulate<F: Flag>(
    type_size: usize,
    f1_vec: &[F],
    f2_vec: &[F],
    g_vec: &[F],
) -> Vec<CsMat<u64>> {
    let k1 = f1_vec[0].size();
    let k2 = f2_vec[0].size();
    let n = g_vec[0].size();
    assert!(type_size <= k1 && k1 <= n);
    assert_eq!(n, k1 + k2 - type_size);
    //
    let shape = (f1_vec.len(), f2_vec.len());
    let storage = if shape.0 < shape.1 { CSR } else { CSC };
    //
    let mut res: Vec<CsMat<u64>> = Vec::with_capacity(g_vec.len());

    for g in g_vec {
        let mut tri_mat = TriMatI::new(shape); //with capacity ?
        let mut iter = Split::with_fixed_part(n, k1, type_size);
        while let Some((subset1, subset2)) = iter.next() {
            let f1 = induce_and_reduce(type_size, g, subset1);
            let f2 = induce_and_reduce(type_size, g, subset2);
            if let Ok(i1) = f1_vec.binary_search(&f1) {
                if let Ok(i2) = f2_vec.binary_search(&f2) {
                    tri_mat.add_triplet(i1, i2, 1);
                }
            }
        }
        let mat = if storage == CSR {
            tri_mat.to_csr()
        } else {
            tri_mat.to_csc()
        };
        res.push(mat)
    }
    res
}

/// Returns a vector `t` where `t[i]` is the index of the unlabeled version
/// of `input_vec[i]`.
///
/// The unlabeling is performed while keeping the image of `eta` as the new type.
pub fn unlabeling_tabulate<F: Flag>(
    eta: &[usize],
    input_vec: &[F],
    output_vec: &[F],
) -> Vec<usize> {
    let type_size = eta.len();
    let mut res = Vec::new();
    for flag in input_vec {
        let unlabeled = canonical_form_typed(&flag.select_type(eta), type_size);
        // FIXME
        if type_size == 0 {
            assert_eq!(&unlabeled, &canonical_form_typed(&unlabeled, 0));
            assert_eq!(&unlabeled, &canonical_form(&unlabeled));
        }
        res.push(
            output_vec.binary_search(&unlabeled).unwrap_or_else(|_| {
                panic!("Flag not found (type {}): {:?}", type_size, &unlabeled)
            }),
        )
    }
    res
}

/// Returns a vector `t` where `t[i]` is the number of ways to relabel
/// `input_vec[i]` while keeping the same (labeled) flag.
///
/// The relabeling is performed while fixing the image of `eta` as the new type.
pub fn unlabeling_count_tabulate<F: Flag>(
    eta: &[usize],
    type_size: usize,
    input_vec: &[F],
) -> Vec<u64> {
    assert!(!input_vec.is_empty());
    let mut res: Vec<u64> = vec![0; input_vec.len()];
    for (i, flag) in input_vec.iter().enumerate() {
        let flag2 = flag.select_type(eta);
        let mut iter = Injection::with_fixed_part(flag2.size(), type_size, eta.len());
        while let Some(phi) = iter.next() {
            let f = canonical_form_typed(&flag2.select_type(phi), type_size);
            if &f == flag {
                res[i] += 1;
            }
        }
    }
    res
}

///Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::*;

    fn count_subflags_ext<F: Flag>(sigma: usize, h: &F, g: &F) -> u64 {
        count_subflags(sigma, &canonical_form_typed(h, sigma), g)
    }

    #[test]
    fn unit_count_subflags() {
        let cherry = Graph::new(3, &[(0, 1), (1, 2)]);
        let c4 = Graph::cycle(4);
        assert_eq!(4, count_subflags_ext(0, &cherry, &c4));
        assert_eq!(2, count_subflags_ext(1, &cherry, &c4));

        let c5 = Graph::cycle(5);
        assert_eq!(
            12,
            count_subflags(0, &canonical_form(&c5), &canonical_form(&Graph::petersen()))
        );
    }

    #[test]
    fn count_subflags_typed() {
        let p3 = Graph::new(3, &[(0, 1), (1, 2)]);
        let p4 = Graph::new(4, &[(0, 1), (1, 2), (2, 3)]);
        let g = Graph::new(6, &[(2, 0), (0, 5), (2, 5), (0, 4), (4, 3)]);
        let g2 = Graph::new(
            7,
            &[
                (0, 1),
                (0, 5),
                (1, 2),
                (1, 4),
                (2, 3),
                (2, 4),
                (5, 4),
                (2, 5),
            ],
        );
        let g3 = Graph::new(5, &[(0, 1), (1, 2), (1, 4), (2, 4), (2, 3)]);

        assert_eq!(1, count_subflags_ext(1, &p3, &g));
        assert_eq!(1, count_subflags_ext(3, &p4, &g2));
        assert_eq!(1, count_subflags_ext(2, &g3, &g2));
    }

    #[test]
    fn unit_count_split() {
        let cherry = canonical_form(&Graph::new(3, &[(0, 1), (1, 2)]));
        let edge = canonical_form(&Graph::new(2, &[(0, 1)]));
        let g = canonical_form(&Graph::new(
            5,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 2)],
        ));
        assert_eq!(4, count_split(0, &cherry, &edge, &g));
    }

    fn csm_get(mat: &CsMat<u64>, i: usize, j: usize) -> u64 {
        assert!(i < mat.rows() && j < mat.cols());
        *mat.get(i, j).unwrap_or(&0)
    }

    fn count_subflag_consistency<F: Flag>(k: usize, n: usize) {
        let a = F::generate(k);
        let b = F::generate(n);
        let tab = count_subflag_tabulate(0, &a, &b);
        for i in 0..a.len() {
            for j in 0..b.len() {
                assert_eq!(count_subflags(0, &a[i], &b[j]), csm_get(&tab, j, i))
            }
        }
    }

    fn count_split_consistency<F: Flag>(k: usize, n: usize) {
        let a = F::generate(k);
        let b = F::generate(n - k);
        let c = F::generate(n);
        let tab = count_split_tabulate(0, &a, &b, &c);
        for i in 0..a.len() {
            for j in 0..b.len() {
                for m in 0..c.len() {
                    assert_eq!(count_split(0, &a[i], &b[j], &c[m]), csm_get(&tab[m], i, j))
                }
            }
        }
    }

    #[test]
    fn unit_count_subflag_tabulate() {
        count_subflag_consistency::<Graph>(3, 5);
    }

    #[test]
    fn unit_count_split_tabulate() {
        count_split_consistency::<Graph>(2, 5);
    }
}
