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
    #[error("Error while parsing matrix coefficient: {0}")]
    ParseError(String),
    #[error("{0}")]
    Io(#[from] io::Error),
    #[error("Solver returned with code {0}")]
    SdpNotSolved(i32),
}

use Error::*;

impl From<ParseIntError> for Error {
    fn from(e: ParseIntError) -> Self {
        ParseError(format!("{}", e))
    }
}

impl From<ParseFloatError> for Error {
    fn from(e: ParseFloatError) -> Self {
        ParseError(format!("{}", e))
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
