use flag_algebra::flags::{OrientedGraph, TriangleFree};
use flag_algebra::*;

type F = SubClass<OrientedGraph, TriangleFree>;
type N = f64;
type V = QFlag<N, F>;

// Parameters
const FLAG_SIZE: usize = 5; // Size of the flags used // can be pushed to 6
const C: f64 = 0.331; // Constant for which we prove the result
const C1: f64 = C; // Constant for which we know it holds

pub fn main() {
    init_default_log();

    // Work basis.
    let b = Basis::new(FLAG_SIZE);

    // 1. Outderee is c.
    let b21 = Basis::new(2).with_type(Type::from_flag(&OrientedGraph::new(1, [])));
    let out_edge = b21.flag(&OrientedGraph::new(2, [(0, 1)]).into());

    let outdegree_is_c = out_edge.at_least(C).multiply_and_unlabel(b).equality();

    // 2. Density of forks is at least 3(3c-1)^2/a

    // Constant for the Chudnovsky-Seymour-Sullivan conjecture
    let a = 0.88;

    // fork
    let fork = flag(&OrientedGraph::new(3, [(0, 1), (0, 2)]).into());

    let fork_ineq = {
        let x = 3. * C - 1.;
        let y = (3. * x * x) / a;
        fork.at_least(y)
    };

    // 3. [|f(sigma)F|] >= 0 for every sigma-source
    // Determine if the flag in input is a sigma-source
    fn is_a_sigma_source(flag: &F, type_size: usize) -> bool {
        for v in type_size..flag.size() {
            // for v not in the type
            if flag.content.out_nbrs(v).iter().any(|&u| u < type_size) {
                return false;
            }
        }
        true
    }

    // Return the sum of sigma sources of the basis in input
    fn sum_of_sigma_sources(basis: Basis<F>) -> V {
        assert_eq!(basis.size, basis.t.size + 1);
        basis.qflag_from_indicator(is_a_sigma_source)
    }
    // Recover the type of the basis in input (as a graph)
    fn get_graph_type(basis: Basis<F>) -> F {
        let graphs: Vec<F> = Basis::new(basis.t.size).get();
        graphs[basis.t.id].clone()
    }

    // Return the vector corresponding to f0
    fn f0(basis: Basis<F>) -> V {
        assert_eq!(basis.size, basis.t.size + 1);
        basis.flag(&get_graph_type(basis).content.add_sink().into())
    }

    // build f(sigma) >= 0
    // where f(sigma) = sum of sigma-sources + (c1-1)F0 - c
    fn f_inequality(basis: Basis<F>) -> Ineq<N, F> {
        let x_f0 = f0(basis) * (C1 - 1.);
        let c_one = basis.one() * C;
        let sum = sum_of_sigma_sources(basis);
        (sum + x_f0 - c_one).at_least(0.)
    }

    fn has_dominated_vertex(g: &F) -> bool {
        for u in 0..g.size() {
            if g.content.in_nbrs(u).len() == g.size() - 1 {
                return true;
            }
        }
        false
    }

    // Build the inequalities of the type f(sigma) >= 0
    //for the types sigma with a dominated vertex
    let mut f_rooted_ineqs = Vec::new();
    for n in 2..b.size {
        for (id, type_flag) in Basis::new(n).get().iter().enumerate() {
            if has_dominated_vertex(type_flag) {
                let basis = Basis::new(n + 1).with_type(Type::new(n, id));
                f_rooted_ineqs.push(f_inequality(basis));
            }
        }
    }
    // Build the inequalities of the type [|f(sigma)*F|] >= 0
    let f_inequalities: Vec<_> = f_rooted_ineqs
        .iter()
        .map(|ineq| ineq.clone().multiply_and_unlabel(b))
        .collect();

    let mut ineqs = f_inequalities;
    ineqs.push(outdegree_is_c);
    ineqs.push(fork_ineq.multiply_and_unlabel(b));
    ineqs.push(total_sum_is_one(b));
    ineqs.push(flags_are_nonnegative(b));

    // Definition of the optimization problem.
    let pb = Problem::<f64, _> {
        // Constraints
        ineqs,
        // Use all relevant Cauchy-Schwarz inequalities.
        cs: b.all_cs(),
        // Minimize density of triangle plus density of independent of size 3.
        obj: b.one(),
    };

    // Write the correspondind SDP program in "cacceta-haggkvist.sdpa".
    // This program can then be solved by CSDP.
    // Infeasible dual means that caccetta-haggkvist is true for d+ = C*n.
    pb.write_sdpa("caccetta-haggkvist").unwrap();
}
