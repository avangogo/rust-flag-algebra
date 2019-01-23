//! Print expression of computations in the flag algebra.

extern crate num;

use std::fmt::*;

/// Expressions that represent a computation in flag algebras.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Expr {
    Add(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>),
    Unlab(Box<Expr>),
    Zero,
    One,
    Num(String),
    Flag(String),
    Var(usize),
}

impl Expr {
    pub fn add(a: Self, b: Self) -> Self {
        Add(Box::new(a), Box::new(b))
    }
    pub fn mul(a: Self, b: Self) -> Self {
        Mul(Box::new(a), Box::new(b))
    }
    pub fn neg(self) -> Self {
        Neg(Box::new(self))
    }
    pub fn sub(a: Self, b: Self) -> Self {
        Self::add(a, b.neg())
    }
    pub fn unlab(self) -> Self {
        Unlab(Box::new(self))
    }
    pub fn num<N>(n: &N) -> Self
    where
        N: num::Num + Display,
    {
        if n == &N::zero() {
            Zero
        } else if n == &N::one() {
            One
        } else {
            Num(format!("{}", n))
        }
    }
}

use self::Expr::*;

impl Display for Expr {
    fn fmt(&self, f: &mut Formatter) -> Result {
        match self.simplify0() {
            Add(a, b) => {
                if let Neg(b1) = *b {
                    write!(f, "{} - {}", a, b1)
                } else {
                    write!(f, "{} + {}", a, b)
                }
            }
            Mul(a, b) => write!(f, "{}*{}", a, b),
            Neg(a) => write!(f, "-{}", a),
            Unlab(a) => write!(f, "[|{}|]", a),
            Zero => write!(f, "0"),
            One => write!(f, "1"),
            Num(s) => write!(f, "{}", s),
            Flag(s) => write!(f, "{}", s),
            Var(_) => write!(f, "x"),
        }
    }
}

impl Expr {
    fn simplify0(&self) -> Self {
        match self {
            Add(a0, b0) => {
                let a = a0.simplify0();
                let b = b0.simplify0();
                if a == Zero {
                    b
                } else if b == Zero {
                    a
                } else {
                    Expr::add(a, b)
                }
            }
            Mul(a0, b0) => {
                let a = a0.simplify0();
                let b = b0.simplify0();
                if a == One {
                    b
                } else if b == One {
                    a
                } else if a == Zero || b == Zero {
                    Zero
                } else {
                    Expr::mul(a, b)
                }
            }
            Neg(a0) => {
                let a = a0.simplify0();
                if a == Zero {
                    Zero
                } else {
                    Expr::neg(a)
                }
            }
            Unlab(a0) => {
                let a = a0.simplify0();
                if a == Zero {
                    Zero
                } else {
                    Expr::unlab(a)
                }
            }
            a => a.clone(),
        }
    }
}
