# flag-algebra

An implementation of
[flag algebras](http://people.cs.uchicago.edu/~razborov/files/flag.pdf).

Flag algebras is a framework used to produce computer-assisted proofs of some inequalities in combinatorics, relying on Semidefinite programming.

## Example

```rust
// Proving that in any graph, at least 1/4 of the triples are triangles or independent sets.
extern crate flag_algebra;

use flag_algebra::*;
use flag_algebra::flags::Graph;
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
   // This program can then be solved by CSDP. The answer would be 0.25.
   pb.write_sdpa("goodman").unwrap();
}
```


License: GPL-3.0
