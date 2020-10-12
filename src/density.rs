//! Computing coefficients of a flag algebra operator.

use crate::flag::Flag;
use crate::iterators::*;
use sprs::{CsMat, TriMatI, CSC, CSR};

fn induce_and_reduce<F: Flag>(type_size: usize, f: &F, subset: &[usize]) -> F {
    f.induce(subset).canonical_typed(type_size)
}

/// Returns the number of induced subflags of `g` isomorphic to `f`,
/// considered with type of size `type_size`.
pub fn count_subflags<F: Flag>(type_size: usize, f: &F, g: &F) -> u32 {
    // f must be in normal form
    assert_eq!(*f, f.canonical_typed(type_size));
    let k = f.size();
    let n = g.size();
    assert!(type_size <= k && k <= n);
    let mut count = 0;
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
pub fn count_split<F: Flag>(type_size: usize, f1: &F, f2: &F, g: &F) -> u32 {
    assert_eq!(*f1, f1.canonical_typed(type_size));
    assert_eq!(*f2, f2.canonical_typed(type_size));
    let k1 = f1.size();
    let k2 = f2.size();
    let n = g.size();
    assert!(type_size <= k1 && k1 <= n);
    assert_eq!(n, k1 + k2 - type_size);
    let mut count = 0;
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
pub fn count_subflag_tabulate<F: Flag>(type_size: usize, f_vec: &[F], g_vec: &[F]) -> CsMat<u32> {
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
            match f_vec.binary_search(&f0) {
                Ok(i) => column[i] += 1,
                Err(_) => {
                    if F::HEREDITARY {
                        panic!("Flag not found: {:?}", &f0)
                    }
                }
            };
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
) -> Vec<CsMat<u32>> {
    let k1 = f1_vec[0].size();
    let k2 = f2_vec[0].size();
    let n = g_vec[0].size();
    assert!(type_size <= k1 && k1 <= n);
    assert_eq!(n, k1 + k2 - type_size);
    //
    let shape = (f1_vec.len(), f2_vec.len());
    let storage = if shape.0 < shape.1 { CSR } else { CSC };
    //
    let mut res: Vec<CsMat<u32>> = Vec::with_capacity(g_vec.len());

    for g in g_vec {
        let mut tri_mat = TriMatI::new(shape); //with capacity ?
        let mut iter = Split::with_fixed_part(n, k1, type_size);
        while let Some((subset1, subset2)) = iter.next() {
            let f1 = induce_and_reduce(type_size, g, subset1);
            let f2 = induce_and_reduce(type_size, g, subset2);
            if let (Ok(i1), Ok(i2)) = (f1_vec.binary_search(&f1), f2_vec.binary_search(&f2)) {
                tri_mat.add_triplet(i1, i2, 1);
            } else if F::HEREDITARY {
                panic!("Flag not found")
            };
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
        let unlabeled = flag.select_type(eta).canonical_typed(type_size);
        // FIXME
        if type_size == 0 {
            debug_assert_eq!(unlabeled, unlabeled.canonical_typed(0));
            debug_assert_eq!(unlabeled, unlabeled.canonical());
        }
        res.push(
            output_vec.binary_search(&unlabeled).unwrap_or_else(|_| {
                panic!("Flag not found (type {}): {:?}", type_size, &unlabeled)
            }),
        )
    }
    res
}

/*
/// Returns a vector `t` where `t[i]` is the number of ways to relabel
/// `input_vec[i]` while keeping the same (labeled) flag.
///
/// The relabeling is performed while fixing the image of `eta` as the new type.
pub fn old_unlabeling_count_tabulate<F: Flag>(
    eta: &[usize],
    type_size: usize,
    input_vec: &[F],
) -> Vec<u32> {
    assert!(!input_vec.is_empty());
    let mut res = vec![0; input_vec.len()];
    for (i, flag) in input_vec.iter().enumerate() {
        let flag2 = flag.select_type(eta).canonical_typed(type_size);
        let mut iter = Injection::with_fixed_part(flag2.size(), type_size, eta.len());
        while let Some(phi) = iter.next() {
            let f = flag2.select_type(phi).canonical_typed(type_size);
            if f == flag2 {
                res[i] += 1;
            }
        }
    }
    res
}
*/

pub fn unlabeling_count_tabulate<F: Flag>(
    eta: &[usize],
    type_size: usize,
    input_vec: &[F],
) -> Vec<u32> {
    // Can be further optimized
    assert!(!input_vec.is_empty());
    let mut res = vec![0; input_vec.len()];
    for (i, flag) in input_vec.iter().enumerate() {
        let flag2 = flag.select_type(eta).canonical_typed(eta.len());
        let aut_typed = flag.stabilizer(type_size).count() as u32;
        let aut = flag2.stabilizer(eta.len()).count() as u32;
        assert_eq!(aut % aut_typed, 0);
        res[i] = aut / aut_typed;
    }
    res
}

/// Classify the typed flags in `input_vec` with type size `type_size` according to their
/// modulo the symmetries of the type.
///
/// The output is a vector `class` where `class[i] == class[j]` if and only if  the rooted
/// flags `input_vec[i]` and `input_vec[j]` are equal modulo permutations of the type.
/// The values in `class` are integers form 0 to "number of classes" - 1 assigned in increasing
/// order.
pub fn invariant_classes<F: Flag>(eta: &[usize], type_size: usize, input_vec: &[F]) -> Vec<usize> {
    // This can be simplified by computing normal form from a partition corresponding to
    // {type, rest} if such a feature is exposed in the canonize crate
    assert!(!input_vec.is_empty());
    let n = input_vec[0].size();
    assert!(type_size <= n);
    if !eta.is_empty() {
        unimplemented!()
    };
    let type_set: Vec<_> = (0..type_size).collect();
    let type_flag = input_vec[0].induce(&type_set).canonical_typed(eta.len());
    // Precompute the symmetry group of the type
    let type_automorphisms: Vec<_> = type_flag
        .stabilizer(eta.len())
        .map(|mut f| {
            f.extend(type_size..n);
            f
        })
        .collect();
    // class[i] always contain Some(c) if the class of input_vec[i] is found and is c
    let mut class = vec![None; input_vec.len()];
    let mut n_classes = 0;
    for (i, flag) in input_vec.iter().enumerate() {
        if class[i] == None {
            let new_class = n_classes;
            n_classes += 1;
            for phi in &type_automorphisms {
                let flag2 = flag.apply_morphism(phi).canonical_typed(type_size);
                let id = input_vec.binary_search(&flag2).expect("Flag not found");
                match class[id] {
                    Some(class_id) => assert_eq!(class_id, new_class),
                    None => class[id] = Some(new_class),
                }
            }
        }
    }
    class.into_iter().map(|c| c.unwrap()).collect()
}

/// Matrices to separate the invariant and anti-invariant parts of a quantum flag.
pub fn class_matrices(class: &[usize]) -> (CsMat<i64>, CsMat<i64>) {
    assert!(!class.is_empty());
    let n_classes = *class.iter().max().unwrap() + 1;
    let n_flags = class.len();
    let mut invariant = TriMatI::new((n_flags, n_classes));
    let mut antiinvariant = TriMatI::new((n_flags, n_flags - n_classes));
    // first_of_class contains the first element of the class i, if met.
    // If class is valid, the classes are encounter in the order of their id
    let mut first_of_class = Vec::new();
    for (i, &c) in class.iter().enumerate() {
        // build the invariant matrix
        invariant.add_triplet(i, c, 1);
        // build the antiinvariant matrix
        if c < first_of_class.len() {
            let col = i - first_of_class.len();
            antiinvariant.add_triplet(i, col, 1);
            antiinvariant.add_triplet(first_of_class[c], col, -1);
        } else {
            assert_eq!(c, first_of_class.len());
            first_of_class.push(i)
        }
    }
    (invariant.to_csc(), antiinvariant.to_csc())
}

///Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::*;
    use canonical_form::Canonize;

    fn count_subflags_ext<F: Flag>(sigma: usize, h: &F, g: &F) -> u32 {
        count_subflags(sigma, &h.canonical_typed(sigma), g)
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
            count_subflags(0, &c5.canonical(), &Graph::petersen().canonical())
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
        let cherry = Graph::new(3, &[(0, 1), (1, 2)]).canonical();
        let edge = Graph::new(2, &[(0, 1)]).canonical();
        let g = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 2)]).canonical();
        assert_eq!(4, count_split(0, &cherry, &edge, &g));
    }

    fn csm_get<N: num::Zero + Copy>(mat: &CsMat<N>, i: usize, j: usize) -> N {
        assert!(i < mat.rows() && j < mat.cols());
        *mat.get(i, j).unwrap_or(&N::zero())
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
        let left = F::generate(k);
        let right = F::generate(n - k);
        let prod = F::generate(n);
        let tab = count_split_tabulate(0, &left, &right, &prod);
        for i in 0..left.len() {
            for j in 0..right.len() {
                for m in 0..prod.len() {
                    assert_eq!(
                        count_split(0, &left[i], &right[j], &prod[m]),
                        csm_get(&tab[m], i, j)
                    )
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
