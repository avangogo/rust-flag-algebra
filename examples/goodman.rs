extern crate flag_algebra;

use flag_algebra::*;
use flags::Graph;
use operator::Basis;
use sdp::Problem;

pub fn main() {
    // Work on the graphs of size 3.
    let basis = Basis::new(3);

    // Define useful flags.
    let k3 = flag(&Graph::new(3, &[(0, 1), (1, 2), (2, 0)])); // Triangle
    let e3 = flag(&Graph::new(3, &[])); // Independent set of size 3

    // Definition of the optimization problem.
    let pb = Problem::<i64, _> {
        // Constraints
        ineqs: vec![total_sum_is_one(basis), flags_are_nonnegative(basis)],
        // Use all relevant Cauchy-Schwarz inequalities.
        cs: basis.all_cs(),
        // Minimize density of triangle plus density of independent of size 3.
        obj: k3 + e3,
    };

    // Write the correspondind SDP program in "goodman.sdpa".
    // This program can then be solved by CSDP.
    pb.write_sdpa("goodman").unwrap();
}
