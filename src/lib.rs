//!A generic implementation of
//![flag algebras](http://people.cs.uchicago.edu/~razborov/files/flag.pdf).
//!
//! Flag algebras is a framework used to produce computer-assisted proofs of some inequalities in combinatorics, relying on Semi-Definite Programming.
//!
//!# Example
//!
//!```rust,no_run
//! // Proving that in any graph, at least 1/4 of the triples
//! // are triangles or independent sets.
//!extern crate flag_algebra;
//!
//!use flag_algebra::*;
//!use flag_algebra::flags::Graph;
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
//!# Features
//! This library can currently do the following.
//! * Generate list of flags from scratch.
//! * Generate flag algebra operators and memoize them in files.
//! * Compute in the flag algebra (multiplication, unlabeling) and add user-defined vectors.
//! * Define, manipulate or amplify flag inequalities (for instance by multiplying an inequality by all flags).
//! * Write problem in .spda format or directly run the CSDP solver.
//! * Automatically eliminate unnecessary constraints (in a naive way).
//! * It is generic:
//! defining new specific class/subclass of flags boils down to implementing a Rust Trait.
//! * Output flags, problems or certificates as html pages
//! in (hopefully) human-readable format (provided that it has a reasonnable size).
//!
//!# Supported flags
//! This library is generic.
//! To use a kind combinatorial objects as flags (e.g. graphs), it suffices to
//! implement the [Flag](trait.Flag.html) trait for the corresponding Rust datatype.
//!
//! Currently, [Flag](trait.Flag.html) is implemented for [Graphs](flags/struct.Graph.html),
//! [Digraphs](flags/struct.Digraph.html) and [edge-colored graphs](flags/struct.CGraph.html)
//! with some fixed number of colors.
//!
//! Beside implementing directly [Flag](trait.Flag.html) for your own types, two mechanisms help
//! to define flag classes based on an existing flag class `F`.
//! * The [Colored](flags/struct.Colored.html) structure for defining vertex-colored flags.
//! If `N` is an integer identifier, `Colored<F, N>` is the type for flags of type `F`
//! where the vertices are further colored in `N` different colors.
//! `Colored<F, N>` automatically implement `Flag` when `F` does.
//! * The [Subclass](struct.SubClass.html) structure and
//! the [SubFlag](trait.SubFlag.html) for classes that are subsets
//! of already defined classes.
//! This is usefull for instance for computing in triangle-free graphs flag algebra
//! without considering other graphs.

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
pub use crate::operator::{Basis, Savable, Type};

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
