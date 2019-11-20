extern crate flag_algebra;

use flag_algebra::*;
use flags::Graph;
use sdp::Problem;

fn is_hamiltonian(g: &Graph) -> bool {
    let mut path = Vec::new();
    fn explore(u: usize, path: &mut Vec<usize>, g: &Graph) -> bool {
        if path.len() == g.size() {
            return u == 0;
        }
        for v in g.nbrs(u) {
            if !path.contains(&v) {
                path.push(v);
                if explore(v, path, g) {
                    return true;
                };
                path.pop();
            }
        }
        false
    };
    explore(0, &mut path, &g)
}

pub fn main() {
    let b = Basis::<Graph>::new(5);
    let flags_with_c5: QFlag<i64, _> =
        b.from_indicator(|g, _| is_hamiltonian(g) );
    let edge = flag(&Graph::new(2, &[]));
    let pb = Problem {
        ineqs: vec![
            total_sum_is_one(b),
            flags_are_nonnegative(b),
            flags_with_c5.at_most(0),
        ],
        cs: b.all_cs(),
        obj: edge.expand(b),
    };

    pb.write_sdpa("c5-free").unwrap();
}
