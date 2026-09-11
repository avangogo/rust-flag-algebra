use std::fmt::Display;
use std::io;
use std::num::{ParseFloatError, ParseIntError};
use std::result::Result;
use std::str::FromStr;
use thiserror::Error;

// A line in a .sdpa format
#[derive(Debug, Clone, Copy)]
pub struct SdpaCoeff {
    pub mat: usize,
    pub block: usize,
    pub i: usize,
    pub j: usize,
    pub val: f64,
}

#[derive(Error, Debug)]
pub enum Error {
    #[allow(clippy::enum_variant_names)]
    #[error("Error while parsing matrix coefficient: {0}")]
    ParseError(String),
    #[error("{0}")]
    Io(#[from] io::Error),
    #[error("Solver did not solve: {0:?}")]
    NotSolved(Outcome),
}

use Error::*;

/// What the solver concluded.
///
/// The verdicts are kept apart because they are different answers: a program
/// whose dual is infeasible is decided, one the solver gave up on is not.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// Solved: `primal` and `dual` bracket the optimum.
    Solved { primal: f64, dual: f64 },
    /// Solved, but not to full accuracy.
    Inaccurate { primal: f64, dual: f64 },
    /// The primal has no feasible point.
    PrimalInfeasible,
    /// The dual has no feasible point.
    DualInfeasible,
    /// The solver stopped without deciding.
    Inconclusive { code: i32, message: String },
}

impl Outcome {
    /// The objective value, when the solver reached one.
    pub fn value(&self) -> Option<f64> {
        match *self {
            Solved { primal, .. } | Inaccurate { primal, .. } => Some(primal),
            _ => None,
        }
    }

    /// The same verdict with the objective values divided by `scale`.
    pub fn scaled(self, scale: f64) -> Self {
        match self {
            Solved { primal, dual } => Solved {
                primal: primal / scale,
                dual: dual / scale,
            },
            Inaccurate { primal, dual } => Inaccurate {
                primal: primal / scale,
                dual: dual / scale,
            },
            other => other,
        }
    }

    /// Read the verdict from csdp's exit code and the lines it ended with.
    pub(crate) fn from_csdp(code: i32, lines: &[String]) -> Result<Self, Error> {
        let objective = |name: &str| -> Result<f64, Error> {
            let prefix = format!("{name} objective value:");
            let line = lines
                .iter()
                .find_map(|line| line.strip_prefix(&prefix))
                .ok_or_else(|| ParseError(format!("csdp printed no {prefix}")))?;
            line.trim().parse().map_err(Error::from)
        };
        Ok(match code {
            0 => Solved {
                primal: objective("Primal")?,
                dual: objective("Dual")?,
            },
            3 => Inaccurate {
                primal: objective("Primal")?,
                dual: objective("Dual")?,
            },
            1 => PrimalInfeasible,
            2 => DualInfeasible,
            _ => Inconclusive {
                code,
                message: lines.first().cloned().unwrap_or_default(),
            },
        })
    }
}

use Outcome::*;

impl From<ParseIntError> for Error {
    fn from(e: ParseIntError) -> Self {
        ParseError(format!("{e}"))
    }
}

impl From<ParseFloatError> for Error {
    fn from(e: ParseFloatError) -> Self {
        ParseError(format!("{e}"))
    }
}

impl FromStr for SdpaCoeff {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut iter = s.split_whitespace();
        let mut next = || {
            iter.next()
                .ok_or_else(|| ParseError("Less than 5 elements".into()))
        };
        let result = SdpaCoeff {
            mat: next()?.parse()?,
            block: next()?.parse()?,
            i: next()?.parse()?,
            j: next()?.parse()?,
            val: next()?.parse()?,
        };
        if iter.next().is_some() {
            return Err(ParseError("Less than 5 elements".into()));
        };
        Ok(result)
    }
}

impl Display for SdpaCoeff {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {} {} {} {}",
            self.mat, self.block, self.i, self.j, self.val
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn outcome(code: i32, lines: &[&str]) -> Outcome {
        let lines: Vec<String> = lines.iter().map(|l| (*l).to_owned()).collect();
        Outcome::from_csdp(code, &lines).unwrap()
    }

    /// The lines are csdp 6.2.0's, verbatim.
    #[test]
    fn csdp_verdicts() {
        assert_eq!(
            outcome(
                0,
                &[
                    "Success: SDP solved",
                    "Primal objective value: 4.9999999e+00 ",
                    "Dual objective value: 5.0000000e+00 ",
                ]
            ),
            Solved {
                primal: 4.9999999,
                dual: 5.0
            }
        );
        assert_eq!(
            outcome(1, &["Success: SDP is primal infeasible"]),
            PrimalInfeasible
        );
        assert_eq!(
            outcome(2, &["Success: SDP is dual infeasible"]),
            DualInfeasible
        );
        assert_eq!(
            outcome(
                3,
                &[
                    "Partial Success: SDP solved with reduced accuracy",
                    "Primal objective value: 1.0e+00 ",
                    "Dual objective value: 2.0e+00 ",
                ]
            ),
            Inaccurate {
                primal: 1.0,
                dual: 2.0
            }
        );
        // Proved and undecided used to be the same `Err`.
        assert_eq!(
            outcome(7, &["Failure: Lack of progress"]),
            Inconclusive {
                code: 7,
                message: "Failure: Lack of progress".to_owned()
            }
        );
        assert_eq!(outcome(7, &["Failure: Lack of progress"]).value(), None);
    }
}
