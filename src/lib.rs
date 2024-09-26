//!A generic implementation of
//![flag algebras](http://people.cs.uchicago.edu/~razborov/files/flag.pdf).
//!
//! Flag algebras is a framework used to produce computer-assisted proofs of some inequalities in combinatorics, relying on Semi-Definite Programming.
//!
//!# Example
//!
//!```rust,no_run
//!// Proving that in any graph, at least 1/4 of the triples
//!// are triangles or independent sets.
//!use flag_algebra::*;
//!use flag_algebra::flags::Graph;
//!
//!// Work on the graphs of size 3.
//!let basis = Basis::new(3);
//!
//!// Define useful flags.
//!let k3 = flag(&Graph::new(3, &[(0, 1), (1, 2), (2, 0)])); // Triangle
//!let e3 = flag(&Graph::new(3, &[])); // Independent set of size 3
//!
//!// Definition of the optimization problem.
//!let pb = Problem::<i64, _> {
//!     // Constraints
//!    ineqs: vec![total_sum_is_one(basis), flags_are_nonnegative(basis)],
//!     // Use all relevant Cauchy-Schwarz inequalities.
//!     cs: basis.all_cs(),
//!     // Minimize density of triangle plus density of independent of size 3.
//!     obj: k3 + e3,
//! };
//!
//! // Write the correspondind SDP program in "goodman.sdpa".
//! // This program can then be solved by CSDP. The answer would be 0.25.
//! pb.write_sdpa("goodman").unwrap();
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
//! implement the [Flag](trait@Flag) trait for the corresponding Rust datatype.
//!
//! Currently, [Flag](trait@Flag) is implemented for [Graphs](struct@flags::Graph),
//! [Oriented graphs](struct@flags::OrientedGraph), [Directed graphs](struct@flags::DirectedGraph)
//! and [edge-colored graphs](struct@flags::CGraph) with some fixed number of colors.
//!
//! Beside implementing directly [Flag] for your own types, two mechanisms help
//! to define flag classes based on an existing flag class `F`.
//! * The [Colored](struct@flags::Colored) structure for defining vertex-colored flags.
//! If `N` is an integer identifier, `Colored<F, N>` is the type for flags of type `F`
//! where the vertices are further colored in `N` different colors.
//! `Colored<F, N>` automatically implement `Flag` when `F` does.
//! * The [`SubClass`] structure and
//! the [`SubFlag`] for classes that are subsets
//! of already defined classes.
//! This is usefull for instance for computing in triangle-free graphs flag algebra
//! without considering other graphs.
//! See the documentation page of [`SubFlag`] for more details.
//!
//!# Expressing elements of a flag algebra
//! See [Type], [Basis] and [`QFlag`].
//!
//! The `Type<F:Flag>` structure identifies a
//! "type" σ in the sense of flag algebras (i.e. a completely labeled flag)
//! is represented by an object.
//! The `Basis<F:Flag>` structure corresponds to a couple (n, σ)
//! and identifies the set of σ-flags of size n.
//! The structure `QFlag`

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
    unused_results
)]

mod algebra;
pub use crate::algebra::*;

mod combinatorics;
mod density;
pub mod flags;
mod iterators;
pub mod operator;
pub mod tools;
pub use crate::operator::{Basis, Savable, Type};

mod expr;
pub mod sdp;
pub mod sdpa;
pub use crate::sdp::Problem;
pub use crate::tools::FlagSolver;

mod flag;
pub use crate::flag::*;

#[macro_use]
extern crate serde_derive;

// Feedback information in the library are sent using env_logger
// This logs require to be initialized
fn logger(level: log::LevelFilter) {
    env_logger::builder()
        .format_module_path(false)
        .format_target(false)
        .filter_level(level)
        .init()
}

/// Initialize the logs to be output to the console.
///
/// In order to be recorded, the logs need to be initialized via this function
/// or any other logger
pub fn init_default_log() {
    logger(log::LevelFilter::Info)
}

/// Initialize the logs to be outputted to the console with detailed information.
pub fn init_debug_log() {
    logger(log::LevelFilter::Trace)
}
