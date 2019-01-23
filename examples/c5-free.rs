extern crate flag_algebra;

use flag_algebra::*;
use crate::sdp::*;
use crate::flags::*;
use canonical_form::Canonize;

fn is_hamiltonian(g: &Graph) -> bool {
    let mut path = Vec::new(); 
    fn explore(u: usize, path: &mut Vec<usize>, g: &Graph) -> bool {
        if path.len() == g.size() {
            return u == 0
        }
        for v in g.nbrs(u) {
            if !path.contains(&v) {
                path.push(v);
                if explore(v, path, g) { return true };
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
        b.from_indicator(|g, _| { if is_hamiltonian(g) { 1 } else { 0 }});
    let edge = Basis::new(2).flag(&Graph::new(2, &[]));
    let pb = Problem {
        ineqs: vec!(total_sum_is_one(b),
                    flags_are_nonnegative(b),
                    flags_with_c5.at_most(0)),
        cs: b.all_cs(),
        obj: edge.expand(b),
    };

    pb.write_sdpa("c5-free").unwrap();
}
