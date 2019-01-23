//! Create and manipulate semi-definite problems.

extern crate ndarray;
extern crate num;
extern crate sprs;

use self::num::*;
use crate::algebra::*;
use crate::flag::Flag;
use crate::operator::*;

use std::fmt::Display;
use std::fs::File;
use std::io::*;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// An optimization problem expressed in flags algebra.
#[derive(Debug)]
pub struct Problem<N, F> {
    /// Set of contraint inequalities.
    pub ineqs: Vec<Ineq<N, F>>,
    /// Set of Cauchy-Schwarz inequlities to be used.
    pub cs: Vec<MulAndUnlabeling<F>>,
    /// Vector to be optimized.
    pub obj: QFlag<N, F>,
}

fn write_coeff<N>(
    file: &mut BufWriter<File>,
    mat_num: usize,
    block_num: usize,
    i: usize,
    j: usize,
    value: N,
) -> Result<()>
where
    N: Display + Zero + PartialEq + Copy,
{
    if value != N::zero() {
        writeln!(
            file,
            "{} {} {} {} {}",
            mat_num,
            block_num + 1,
            i + 1,
            j + 1,
            value
        )?;
    }
    Ok(())
}

impl<N, F> Problem<N, F>
where
    N: Display + Zero + Copy + PartialEq,
    F: Flag,
{
    /// Panic if the size of the basis involved are inconsistent.
    pub fn check(&self) {
        let b = self.obj.basis;
        for ineq in self.ineqs.iter() {
            assert_eq!(ineq.meta.basis, b);
        }
        for cs in self.cs.iter() {
            assert_eq!(cs.output_basis(), b);
        }
    }

    fn write_header(&self, file: &mut BufWriter<File>) -> Result<()> {
        writeln!(file, "* sdp problem generated by Rust")?;
        writeln!(file, "* Flag: {}", F::name())?;
        writeln!(file, "* Basis: {:}", self.obj.basis)?;
        writeln!(file, "* {} groups of linear constraints", self.ineqs.len())?;
        writeln!(file, "* {} Cauchy-Schwarz constraints", self.cs.len())?;
        for i in self.ineqs.iter() {
            writeln!(
                file,
                "* {} >= {} ({})",
                i.meta.flag_expr,
                i.meta.bound_expr,
                i.data.len()
            )?;
        }
        for cs in self.cs.iter() {
            writeln!(file, "* {}", cs)?;
        }
        writeln!(file, "*")
    }
    /// Write the semi-definite program in the file `filename` in the sdpa format.
    pub fn write_sdpa(&self, filename: &str) -> Result<()> {
        self.check();
        let mut filename = PathBuf::from(filename);
        let _ = filename.set_extension("sdpa");
        let mut file = BufWriter::new(File::create(&filename)?);
        self.write_header(&mut file)?;
        let cs_mat: Vec<_> = self.cs.iter().map(|x| x.get()).collect();
        // Line 1: Number of constraints = size of the basis
        writeln!(file, "{}", self.obj.data.len())?;
        // Line 2: Number of blocks (one for each constraint)
        writeln!(file, "{}", self.ineqs.len() + self.cs.len())?;
        // Line 3: Sizes of the blocks
        for ineq in self.ineqs.iter() {
            write!(file, "-{} ", ineq.data.len())?;
        }
        for split in cs_mat.iter() {
            write!(file, "{} ", split[0].rows())?;
        }
        writeln!(file)?;
        // Line 4: vector ai
        for v in self.obj.data.iter() {
            write!(file, "{} ", v)?;
        }
        writeln!(file)?;
        // Lines 5+: body
        // Matrix 0: Objective
        for (block_num, ineq) in self.ineqs.iter().enumerate() {
            for (i, ref inequality) in ineq.data.iter().enumerate() {
                write_coeff(&mut file, 0, block_num, i, i, inequality.bound)?;
            }
        }
        writeln!(file)?;
        // Matrices 1+:
        // Inequaltity blocks
        for (block_num, ineq) in self.ineqs.iter().enumerate() {
            for (i, ref ineqdata) in ineq.data.iter().enumerate() {
                for (mat_num, &v) in ineqdata.flag.iter().enumerate() {
                    write_coeff(&mut file, mat_num + 1, block_num, i, i, v)?;
                }
            }
        }
        writeln!(file)?;
        // Cs blocks
        let offset = self.ineqs.len();
        for (block_num, line) in cs_mat.iter().enumerate() {
            for (mat_num, matrix) in line.iter().enumerate() {
                for (&v, (i, j)) in matrix.iter() {
                    if i <= j {
                        write_coeff(&mut file, mat_num + 1, block_num + offset, i, j, v)?;
                    }
                }
            }
        }
        Ok(())
    }
}

// ========== Certificate ==========
#[derive(Debug, Clone)]
pub struct Certificate {
    y: Vec<f64>,
}

impl Certificate {
    pub fn from_file(name: &str) -> Result<Self> {
        let file = File::open(name)?;
        let mut buf = BufReader::new(file).lines();
        let line: String = buf.next().unwrap().unwrap();
        let y: Vec<f64> = line[..]
            .split(' ')
            .filter(|x| x.is_empty())
            .map(|x| {
                println!("{}", x);
                x.parse().unwrap()
            })
            .collect();
        Ok(Certificate { y })
    }
    pub fn to_vec<F>(&self, b: Basis<F>, threshold: f64) -> QFlag<f64, F>
    where
        F: Flag,
    {
        b.from_vec(
            self.y
                .iter()
                .map(|x| if x.abs() < threshold { 0. } else { *x })
                .collect(),
        )
    }
}
