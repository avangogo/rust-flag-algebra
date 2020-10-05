//!An implementation of
//![flag algebras](http://people.cs.uchicago.edu/~razborov/files/flag.pdf).
//!
//! Flag algebras is a framework used to produce computer-assisted proofs of some inequalities in combinatorics, relying on Semidefinite programming.
//!
//!# Example
//!
//!```rust,no_run
//! // Proving that in any graph, at least 1/4 of the triples are triangles or independent sets.
//!extern crate flag_algebra;
//!
//!use flag_algebra::*;
//!use flag_algebra::flags::Graph;
//!use sdp::Problem;
//!
//!pub fn main() {
//!    // Work on the graphs of size 3.
//!    let basis = Basis::new(3);
//!
//!    // Define useful flags.
//!    let k3 = flag(&Graph::new(3, &[(0, 1), (1, 2), (2, 0)])); // Triangle
//!    let e3 = flag(&Graph::new(3, &[])); // Independent set of size 3
//!
//!    // Definition of the optimization problem.
//!    let pb = Problem::<i64, _> {
//!        // Constraints
//!        ineqs: vec![total_sum_is_one(basis), flags_are_nonnegative(basis)],
//!        // Use all relevant Cauchy-Schwarz inequalities.
//!        cs: basis.all_cs(),
//!        // Minimize density of triangle plus density of independent of size 3.
//!        obj: k3 + e3,
//!    };
//!
//!    // Write the correspondind SDP program in "goodman.sdpa".
//!    // This program can then be solved by CSDP. The answer would be 0.25.
//!    pb.write_sdpa("goodman").unwrap();
//!}
//!```
//!

#![warn(
    missing_debug_implementations,
    missing_copy_implementations,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_import_braces,
    unused_qualifications,
    unused_labels,
//    unused_results
)]

mod algebra;
pub use crate::algebra::*;

mod combinatorics;
mod common;
pub mod density;
pub mod draw;
pub mod flags;
mod iterators;
pub mod operator;
pub use crate::operator::{Basis, Type, Savable};

pub mod expr;
mod reduction;
pub mod report;
pub use crate::report::Html;
pub mod sdp;
pub use crate::sdp::Problem;
pub mod sdpa;
pub use crate::reduction::*;

mod flag;
pub use crate::flag::*;

#[macro_use]
extern crate serde_derive;

// Feedback information in the library are sent using simplelog
// This logs require to be initialized
use simplelog::*;

/// Initialize the logs to be outputted to the console.
///
/// In order to be recorded, the logs need to be initialized via this function
/// or any initializer of the simplelog library
pub fn init_default_log() {
    let config = ConfigBuilder::new()
        .set_max_level(LevelFilter::Error)
        .set_target_level(LevelFilter::Off)
        .set_thread_level(LevelFilter::Off)
        .build();
    TermLogger::init(LevelFilter::Info, config, TerminalMode::Mixed).unwrap();
}

/// Initialize the logs to be outputted to the console with detailed information.
///
/// In order to be recorded, the logs need to be initialized via this function
/// or any initializer of the simplelog library
pub fn init_debug_log() {
    let config = ConfigBuilder::new()
        .set_max_level(LevelFilter::Error)
        .set_target_level(LevelFilter::Off)
        .set_thread_level(LevelFilter::Off)
        .build();
    TermLogger::init(LevelFilter::Trace, config, TerminalMode::Mixed).unwrap();
}
