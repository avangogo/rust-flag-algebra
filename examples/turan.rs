extern crate flag_algebra;

use flag_algebra::*;
use crate::sdp::Problem;
use crate::flags::Graph;
use crate::operator::Basis;

pub fn main() {
    // Work on the graphs of size 3.
    let basis = Basis::new(3);

    // Define useful flags.
    let triangle = flag(&Graph::new(3, &[(0,1),(1,2),(2,0)]));
    let edge = flag(&Graph::new(2, &[(0,1)]));

    // Definition of the optimization problem.
    let pb = Problem::<i64, _> {
        // Constraints
        ineqs: vec!(
            total_sum_is_one(basis),
            flags_are_nonnegative(basis),
            triangle.at_most(0)
        ),
        // Use all relevant Cauchy-Schwarz inequalities.
        cs: basis.all_cs(),
        // Minimize minus density of edges expressed in the basis of flags of
        // size 3.
        obj: -edge.expand(basis),
    };

    // Write the correspondind SDP program in "turan.sdpa".
    // This program can then be solved by CSDP.
    pb.write_sdpa("turan").unwrap();
}
