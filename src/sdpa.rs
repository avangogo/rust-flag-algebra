use failure::Fail;
use std::fmt::Display;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::num::{ParseFloatError, ParseIntError};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::result::Result;
use std::str::FromStr;

use log::*;

const CS_COST: f64 = 100.;
const INEQ_COST: f64 = 1.;

pub const CERTIFICATE_FILE: &str = "certificate";
pub const CERTIFICATE_MINIMIZE_FILE: &str = "certificate_minimize";

// A line in a .sdpa format
#[derive(Debug, Clone, Copy)]
pub struct SdpaCoeff {
    pub mat: usize,
    pub block: usize,
    pub i: usize,
    pub j: usize,
    pub val: f64,
}

#[derive(Fail, Debug)]
pub enum Error {
    #[fail(display = "Error while parsing matrix coefficient: {}", _0)]
    ParseError(String),
    #[fail(display = "{}", _0)]
    Io(#[cause] io::Error),
    #[fail(display = "Solver returned with code {}", _0)]
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

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Io(e)
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

// SDPA format for problems
// 1. dimension ( =b.len() )
// 2. n_blocks ( =block_sizes.len() )
// 3. block_sizes of length nblock
// 4. b of length dim
// 5+. list of coefficients
#[derive(Debug, Clone)]
pub struct SdpaProblem {
    block_sizes: Vec<i32>,
    b: Vec<f64>,
    coeffs: Vec<SdpaCoeff>,
}

impl SdpaProblem {
    pub fn write(&self, filename: &str) -> Result<(), io::Error> {
        let mut w = BufWriter::new(File::create(&filename)?);
        writeln!(w, "{}", self.b.len())?;
        writeln!(w, "{}", self.block_sizes.len())?;
        for i in &self.block_sizes {
            write!(w, "{} ", i)?;
        }
        writeln!(w)?;
        for x in &self.b {
            write!(w, "{} ", x)?;
        }
        writeln!(w)?;
        for coeff in &self.coeffs {
            writeln!(w, "{} ", coeff)?;
        }
        Ok(())
    }
    pub fn load(filename: &str) -> Result<Self, Error> {
        let mut filename = PathBuf::from(filename);
        let _ = filename.set_extension("sdpa");
        let buf = BufReader::new(File::open(filename)?);
        let mut lines = buf.lines().map(|line| line.unwrap()).filter(|line| {
            let l = line.trim_start();
            !l.starts_with('*') && !l.is_empty()
        });
        let dim: usize = lines.next().unwrap().parse().unwrap();
        let nblock: usize = lines.next().unwrap().parse()?;
        let block_sizes: Vec<i32> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(block_sizes.len(), nblock);
        let b: Vec<f64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(b.len(), dim);
        let mut coeffs: Vec<SdpaCoeff> = Vec::new();
        for line in lines {
            coeffs.push(line.parse().expect(&line))
        }
        Ok(SdpaProblem {
            block_sizes,
            b,
            coeffs,
        })
    }
    // Write the problem of minimizing the cetificate, under the constrint obj=target_value
    pub fn to_certificate_minimization(mut self, target_value: f64) -> Self {
        self.b.push(target_value);
        let new_mat = self.b.len();
        // The old objective becomes the new matrix
        for coeff in &mut self.coeffs {
            if coeff.mat == 0 {
                coeff.mat = new_mat
            }
        }
        // The new objective is the weight of the trace
        push_identities(&mut self.coeffs, 0, &self.block_sizes, -INEQ_COST, -CS_COST);
        self
    }
}

// SDPA format for certificates (as given by csdp)
// 1. vector y
// 2+. list of coefficients for Z and X
// matrix 1: Z, matrix 2: X
#[derive(Debug, Clone)]
pub struct SdpaCertificate {
    y: Vec<f64>,
    coeffs: Vec<SdpaCoeff>,
}

impl SdpaCertificate {
    pub fn write(&self, filename: &str) -> Result<(), io::Error> {
        let filename = PathBuf::from(filename);
        //let _ = filename.set_extension("cert.sdpa");
        let mut w = BufWriter::new(File::create(&filename)?);
        for v in &self.y {
            write!(w, "{} ", v)?;
        }
        writeln!(w)?;
        for coeff in &self.coeffs {
            writeln!(w, "{} ", coeff)?;
        }
        Ok(())
    }
    pub fn load(filename: &str) -> Result<Self, io::Error> {
        let buf = BufReader::new(File::open(filename)?);
        let mut lines = buf
            .lines()
            .map(|line| line.unwrap())
            .filter(|line| !line.trim_start().is_empty());
        let y: Vec<f64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let mut coeffs: Vec<SdpaCoeff> = Vec::new();
        for line in lines {
            coeffs.push(line.parse().expect(&line))
        }
        Ok(SdpaCertificate { y, coeffs })
    }
    pub fn to_certificate_minimization(mut self, problem: &SdpaProblem) -> Self {
        // Add a dimension to y and set y to 0
        self.y.push(-1.);
        //      self.y.push(0.);
        //     for v in &mut self.y {
        //         *v = 0.
        //     }
        // Discard Z
        //        self.coeffs.retain(|coeff|{ coeff.mat == 2 });
        // Write a new Z as -C from the problem
        for coeff in &problem.coeffs {
            if coeff.mat == 0 {
                let new_coeff = SdpaCoeff {
                    mat: 1,
                    val: -coeff.val,
                    ..*coeff
                };
                self.coeffs.push(new_coeff)
            }
        }
        sum_duplicates(&mut self.coeffs);
        self
    }
    // Discard the dual values as they are now invalid
    pub fn from_certificate_minimization(mut self) -> Self {
        // Remove the extra dimension of y and set y to 0
        let _ = self.y.pop().unwrap();
        self.y.push(0.);
        for v in &mut self.y {
            *v = 0.
        }
        // Discard Z
        self.coeffs.retain(|coeff| coeff.mat == 2);
        self
    }
}

/// Run csdp and parse its output
pub fn csdp(filename: &str, initial_solution: Option<&str>) -> Result<f64, Error> {
    let mut command = Command::new("csdp");
    command.arg(filename).arg(CERTIFICATE_FILE);
    if let Some(sol) = initial_solution {
        command.arg(sol);
    };
    info!("Calling CSDP");
    debug!("command: {:?}", command);
    let mut child = command
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to call csdp");
    let output = BufReader::new(child.stdout.take().unwrap());
    //    let output = child.wait_with_output().unwrap();
    //let output = command.output().expect("Failed to call csdp");
    //    eprint!("{}", String::from_utf8(output.stderr).unwrap());
    use std::time::{Duration, Instant};
    let time_start = Instant::now();
    let mut lines = output.lines();
    let mut stream = false;
    let time_before_stream = Duration::from_secs(2);
    while let Some(line_result) = lines.next() {
        let line = line_result.unwrap();
        if line.starts_with("Iter") {
            if !stream && Instant::now() - time_start > time_before_stream {
                stream = true;
                info!(
                    "csdp is taking more than {}s, start streaming output",
                    time_before_stream.as_secs_f32()
                )
            }
            if stream {
                info!("{}", line)
            } else {
                debug!("{}", line)
            }
        } else {
            let code = child
                .wait()
                .expect("csdp wasn't running")
                .code()
                .expect("No exit code");
            if code == 0 {
                let value: f64 = lines
                    .next()
                    .unwrap()
                    .unwrap()
                    .split_whitespace()
                    .nth(3)
                    .unwrap()
                    .parse()
                    .unwrap();
                info!("{} with primal value {}", line, value);
                return Ok(value);
            } else {
                info!("{}", line);
                return Err(SdpNotSolved(code));
            }
        }
    }
    panic!("CSDP output incorrectly parsed")
}

pub fn csdp_minimize_certificate(
    filename: &str,
    initial_solution: Option<&str>,
) -> Result<f64, Error> {
    let val = csdp(filename, initial_solution)?;
    let problem = SdpaProblem::load(filename)?.to_certificate_minimization(val);
    let cert = SdpaCertificate::load(CERTIFICATE_FILE)?.to_certificate_minimization(&problem);
    let filename_minimize = format!("{}.minimize", filename);
    problem.write(&filename_minimize)?;
    cert.write(CERTIFICATE_MINIMIZE_FILE)?;
    info!("Try to minimize certificate");
    match csdp(&filename_minimize, Some(CERTIFICATE_MINIMIZE_FILE)) {
        Ok(_) => {
            info!("Certificate minimized");
            SdpaCertificate::load(CERTIFICATE_MINIMIZE_FILE)?
                .from_certificate_minimization()
                .write(CERTIFICATE_FILE)?;
        }
        Err(_) => {
            warn!("Cannot minimize certificate");
        }
    }
    Ok(val)
}

impl SdpaCoeff {
    fn indices(&self) -> (usize, usize, usize, usize) {
        (self.mat, self.block, self.i, self.j)
    }
}

fn sum_duplicates(coeffs: &mut Vec<SdpaCoeff>) {
    if coeffs.len() >= 2 {
        coeffs.sort_by_key(SdpaCoeff::indices);
        let mut write = 0;
        for read in 1..coeffs.len() {
            if coeffs[read].indices() == coeffs[write].indices() {
                coeffs[write].val += coeffs[read].val
            } else {
                write += 1;
                coeffs.swap(write, read);
            }
        }
        coeffs.truncate(write + 1)
    }
}

fn push_identities(
    coeffs: &mut Vec<SdpaCoeff>,
    matrix_number: usize,
    block_sizes: &[i32],
    scale_diag: f64,
    scale_nondiag: f64,
) {
    for (block, &blocksize) in block_sizes.iter().enumerate() {
        let val = if blocksize < 0 {
            scale_diag
        } else {
            scale_nondiag
        };
        for i in 0..blocksize.abs() as usize {
            coeffs.push(SdpaCoeff {
                mat: matrix_number,
                block: block + 1,
                i: i + 1,
                j: i + 1,
                val,
            })
        }
    }
}
