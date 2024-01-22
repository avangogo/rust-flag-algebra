//! Expression of computations in the flag algebra for prettyprinting.

use crate::operator::{Basis, Savable, Type};
use std::collections::BTreeMap;
use std::fmt::*;
use std::ops::{Add, Mul, Neg, Sub};
use std::rc::Rc;

/// Expressions that represent a computation in flag algebras.
pub enum Expr<N, F: Flag> {
    Add(RcExpr<N, F>, RcExpr<N, F>),
    Mul(RcExpr<N, F>, RcExpr<N, F>),
    Neg(RcExpr<N, F>),
    Unlab(RcExpr<N, F>),
    Zero,
    One,
    Num(Rc<N>),
    Named(RcExpr<N, F>, Rc<String>, bool),
    Var(usize),
    Flag(usize, Basis<F>),
    FromFunction(Rc<dyn Fn(&F, usize) -> N>, Basis<F>),
    FromIndicator(fn(&F, usize) -> bool, Basis<F>),
    Unknown,
}

// Straightforward trait implementations
// (derive is too conservative when working with Rc or PhantomData,
// follow https://github.com/rust-lang/rust/issues/26925)
impl<N, F: Flag> Clone for Expr<N, F> {
    fn clone(&self) -> Self {
        match self {
            Add(a, b) => Add(a.clone(), b.clone()),
            Mul(a, b) => Mul(a.clone(), b.clone()),
            Neg(a) => Neg(a.clone()),
            Unlab(a) => Unlab(a.clone()),
            Num(a) => Num(a.clone()),
            Var(a) => Var(*a),
            Named(a, b, c) => Named(a.clone(), b.clone(), *c),
            Flag(a, b) => Flag(*a, *b),
            FromFunction(a, b) => FromFunction(a.clone(), *b),
            FromIndicator(a, b) => FromIndicator(*a, *b),
            Unknown => Unknown,
            Zero => Zero,
            One => One,
        }
    }
}

#[derive(Debug, Clone)]
pub enum VarRange<F: Flag> {
    InBasis(Basis<F>),
}

#[derive(Debug, Clone)]
pub struct Names<N, F: Flag> {
    pub flags: BTreeMap<(usize, Basis<F>), String>,
    pub types: BTreeMap<Type<F>, String>,
    pub functions: Vec<(String, QFlag<N, F>)>,
    pub sets: Vec<(String, Basis<F>, Vec<F>)>,
}

impl<N, F: Flag> Default for Names<N, F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<N, F: Flag> Names<N, F> {
    pub fn new() -> Self {
        Self {
            flags: BTreeMap::new(),
            types: BTreeMap::new(),
            functions: Vec::new(),
            sets: Vec::new(),
        }
    }
    pub fn is_empty(&self) -> bool {
        self.flags.is_empty()
            && self.types.is_empty()
            && self.functions.is_empty()
            && self.sets.is_empty()
    }
    fn name_flag(&mut self, i: usize, basis: Basis<F>) -> String
    where
        F: Ord,
    {
        self.flags
            .entry((i, basis))
            .or_insert_with(|| format!("F_{{{}}}^{{{}}}", i, basis.print_concise()))
            .clone()
    }
    fn name_type(&mut self, t: Type<F>) -> String {
        let i = self.types.len();
        self.types
            .entry(t)
            .or_insert_with(|| {
                if i == 0 {
                    "\\sigma".to_string()
                } else {
                    format!("\\sigma_{}", i)
                }
            })
            .clone()
    }
    fn name_set(&mut self, f: fn(&F, usize) -> bool, basis: Basis<F>) -> String
    where
        F: Flag,
    {
        let name = format!("S_{}", self.sets.len() + 1);
        let mut set = basis.get();
        set.retain(|x| f(x, basis.t.size));
        self.sets.push((name.clone(), basis, set));
        name
    }
    fn name_function(&mut self, f: Rc<dyn Fn(&F, usize) -> N>, basis: Basis<F>) -> String
    where
        F: Flag,
    {
        let name = format!("f_{}", self.functions.len() + 1);
        self.functions.push((name.clone(), basis.qflag_from_coeff_rc(f)));
        name
    }
}

use Expr::*;
use VarRange::*;

impl<F: Flag> VarRange<F> {
    fn eval<N>(&self, i: usize) -> Expr<N, F> {
        match self {
            InBasis(basis) => Flag(i, *basis),
        }
    }
    pub(crate) fn latex<N>(&self, names: &mut Names<N, F>) -> String {
        match self {
            InBasis(basis) => format!("\\forall H\\in {},\\quad ", latex_basis(basis, names)),
        }
    }
}

impl<N, F: Flag> Add for Expr<N, F> {
    type Output = Self;

    fn add(self, b: Self) -> Self {
        Add(Rc::new(self), Rc::new(b))
    }
}

impl<N, F: Flag> Neg for Expr<N, F> {
    type Output = Self;

    fn neg(self) -> Self {
        Neg(Rc::new(self))
    }
}

impl<N, F: Flag> Sub for Expr<N, F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl<N, F: Flag> Mul for Expr<N, F> {
    type Output = Self;

    fn mul(self, b: Self) -> Self {
        Mul(Rc::new(self), Rc::new(b))
    }
}

type RcExpr<N, F> = Rc<Expr<N, F>>;

impl<N, F: Flag> Expr<N, F> {
    pub fn unlab(self) -> Self {
        Unlab(Rc::new(self))
    }
    pub fn named(self, name: String) -> Self {
        Named(Rc::new(self), Rc::new(name), false)
    }
    pub fn unknown(name: String) -> Self {
        Unknown.named(name)
    }
    pub fn num(n: &N) -> Self
    where
        N: num::Num + Clone,
    {
        if n == &N::zero() {
            Zero
        } else if n == &N::one() {
            One
        } else {
            Num(Rc::new(n.clone()))
        }
    }
    fn simplify(&self) -> Self
    where
        Expr<N, F>: Clone,
    {
        match self {
            Add(a0, b0) => match (a0.simplify(), b0.simplify()) {
                (Zero, a) | (a, Zero) => a,
                (a, b) => a + b,
            },
            Mul(a0, b0) => match (a0.simplify(), b0.simplify()) {
                (One, a) | (a, One) => a,
                (Zero, _) | (_, Zero) => Zero,
                (a, b) => a * b,
            },
            Neg(a0) => match a0.simplify() {
                Zero => Zero,
                a => -a,
            },
            Unlab(a0) => match a0.simplify() {
                Zero => Zero,
                a => Self::unlab(a),
            },
            a => a.clone(),
        }
    }
    fn is_sum(&self) -> bool {
        matches!(self, Add(_, _))
    }
    pub fn latex(&self, names: &mut Names<N, F>) -> String
    where
        N: Display,
        F: Ord,
    {
        self.simplify().latex0(names)
    }
    fn latex0(&self, names: &mut Names<N, F>) -> String
    where
        N: Display,
        F: Ord,
    {
        match self {
            Add(a, b) => {
                if let Neg(b1) = &**b {
                    format!("{} - {}", a.latex0(names), Paren(b1).latex(names))
                } else {
                    format!("{} + {}", a.latex0(names), b.latex0(names))
                }
            }
            Mul(a, b) => format!("{}\\cdot {}", Paren(a).latex(names), Paren(b).latex(names)),
            Neg(a) => format!("-{}", Paren(a).latex(names)),
            Unlab(a) => format!(
                "\\left[\\!\\!\\left[{}\\right]\\!\\!\\right]",
                a.latex0(names)
            ),
            Zero => "0".into(),
            One => "1".into(),
            Num(s) => format!("{}", s),
            Var(_) => "H".into(),
            Named(e, name, latex) => {
                if *latex {
                    format!("\\textrm{{{}}}", name)
                } else {
                    e.latex0(names)
                }
            }
            Flag(i, basis) => names.name_flag(*i, *basis),
            FromFunction(f, b) => format!(
                "\\sum_{{F\\in{}}} {}(F)F",
                latex_basis(b, names),
                names.name_function(f.clone(), *b)
            ),
            FromIndicator(f, b) => format!(
                "\\sum_{{F\\in {}\\subseteq{}}}F",
                names.name_set(*f, *b),
                latex_basis(b, names)
            ),
            Unknown => "Unknown".into(),
        }
    }
}

fn latex_basis<N, F: Flag>(basis: &Basis<F>, names: &mut Names<N, F>) -> String {
    if basis.t.is_empty() {
        format!("\\mathcal{{F}}_{{{}}}", basis.size)
    } else {
        format!(
            "\\mathcal{{F}}^{{{}}}_{{{}}}",
            names.name_type(basis.t),
            basis.size
        )
    }
}

struct Paren<'a, N, F: Flag>(&'a Expr<N, F>);

impl<'a, N, F: Flag> Display for Paren<'a, N, F>
where
    Expr<N, F>: Display,
{
    fn fmt(&self, f: &mut Formatter) -> Result {
        if self.0.is_sum() {
            write!(f, "({})", self.0)
        } else {
            write!(f, "{}", self.0)
        }
    }
}

impl<'a, N, F> Paren<'a, N, F>
where
    N: Display,
    F: Ord + Flag,
{
    fn latex(&self, names: &mut Names<N, F>) -> String {
        if self.0.is_sum() {
            format!("\\left({}\\right)", self.0.latex0(names))
        } else {
            self.0.latex0(names)
        }
    }
}

impl<N, F: Flag> Display for Expr<N, F>
where
    N: Display,
{
    fn fmt(&self, f: &mut Formatter) -> Result {
        match self.simplify() {
            Add(a, b) => {
                if let Neg(b1) = &*b {
                    write!(f, "{} - {}", a, Paren(b1))
                } else {
                    write!(f, "{} + {}", a, b)
                }
            }
            Mul(a, b) => write!(f, "{}*{}", Paren(&a), Paren(&b)),
            Neg(a) => write!(f, "-{}", Paren(&a)),
            Unlab(a) => write!(f, "[|{}|]", a),
            Zero => write!(f, "0"),
            One => write!(f, "1"),
            Num(s) => write!(f, "{}", s),
            Var(_) => write!(f, "x"),
            Named(_, name, _) => write!(f, "{}", name),
            Flag(i, basis) => write!(f, "flag({}:{})", i, basis.print_concise()),
            FromFunction(_, _) => write!(f, "Σ f(F)F"),
            FromIndicator(_, _) => write!(f, "Σ F"),
            Unknown => write!(f, "unknown"),
        }
    }
}

impl<N, F> Debug for Expr<N, F>
where
    F: Flag + Debug,
    N: Debug,
{
    fn fmt(&self, f: &mut Formatter) -> Result {
        match self {
            Add(a, b) => write!(f, "Add({:?}, {:?})", a, b),
            Mul(a, b) => write!(f, "Mul({:?}, {:?})", a, b),
            Named(a, b, c) => write!(f, "Named({:?}, {:?}, {:?})", a, b, c),
            Flag(a, b) => write!(f, "Flag({:?}, {:?})", a, b),
            Neg(a) => write!(f, "Neg({:?})", a),
            Unlab(a) => write!(f, "Unlab({:?})", a),
            Num(a) => write!(f, "Num({:?})", a),
            Var(a) => write!(f, "Var({:?})", a),
            FromFunction(_, b) => write!(f, "FromFunction(_, {:?})", b),
            FromIndicator(_, b) => write!(f, "FromIndicator(_, {:?})", b),
            Unknown => write!(f, "Unknown"),
            Zero => write!(f, "Zero"),
            One => write!(f, "One"),
        }
    }
}

use crate::Flag;
/// Expression evaluation
use crate::QFlag;
use ndarray::ScalarOperand;
use num::FromPrimitive;

#[derive(Clone, Debug)]
enum Val<N, F: Flag> {
    Num(N),
    QFlag(QFlag<N, F>),
}

impl<N, F> Val<N, F>
where
    N: num::Num + Clone + Neg<Output = N>,
    F: Flag,
{
    fn unwrap_qflag(self) -> QFlag<N, F> {
        if let Self::QFlag(qflag) = self {
            qflag
        } else {
            panic!("QFlag expected")
        }
    }
    fn neg(self) -> Self {
        match self {
            Self::Num(n) => Self::Num(-n),
            Self::QFlag(qflag) => Self::QFlag(-&qflag),
        }
    }
}

impl<N, F> Expr<N, F>
where
    N: num::Num + Neg<Output = N> + Clone + FromPrimitive + ScalarOperand + Display,
    F: Flag,
{
    pub fn eval(&self) -> QFlag<N, F> {
        self.eval0(None).unwrap_qflag()
    }
    pub fn eval_with_context(&self, range: &VarRange<F>, id: usize) -> QFlag<N, F> {
        self.eval0(Some((range, id))).unwrap_qflag()
    }
    fn eval0(&self, context: Option<(&VarRange<F>, usize)>) -> Val<N, F> {
        match self {
            Add(a, b) => match (a.eval0(context), b.eval0(context)) {
                (Val::Num(n1), Val::Num(n2)) => Val::Num(n1 + n2),
                (Val::QFlag(f), Val::QFlag(g)) => Val::QFlag(f + g),
                (Val::QFlag(f), Val::Num(n)) | (Val::Num(n), Val::QFlag(f)) => {
                    assert!(F::HEREDITARY);
                    let one = f.basis.one();
                    Val::QFlag(f + one * n)
                }
            },
            Mul(a, b) => match (a.eval0(context), b.eval0(context)) {
                (Val::Num(n1), Val::Num(n2)) => Val::Num(n1 * n2),
                (Val::QFlag(f), Val::QFlag(g)) => Val::QFlag(f * g),
                (Val::Num(n), Val::QFlag(g)) | (Val::QFlag(g), Val::Num(n)) => Val::QFlag(g * n),
            },
            Neg(e) => e.eval0(context).neg(),
            Unlab(e) => Val::QFlag(e.eval0(context).unwrap_qflag().untype()),
            Num(x) => Val::Num((**x).clone()),
            Var(_) => match context {
                Some((range, id)) => range.eval(id).eval0(None),
                None => panic!("Cannot evaluate variable"),
            },
            Named(e, _, _) => e.eval0(context),
            Flag(i, basis) => Val::QFlag(basis.flag_from_id(*i)),
            FromIndicator(f, basis) => Val::QFlag(basis.qflag_from_indicator(*f)),
            FromFunction(f, basis) => Val::QFlag(basis.qflag_from_coeff_rc(f.clone())),
            Zero => Val::Num(N::zero()),
            One => Val::Num(N::one()),
            Unknown => panic!("Cannot evaluate unknown"),
        }
    }
}
impl<N, F> Expr<N, F>
where
    N: Clone,
    F: Flag,
{
    pub fn substitute_option(&self, range_opt: &Option<VarRange<F>>, id: usize) -> Self {
        match range_opt {
            Some(range) => self.substitute(range, id),
            None => self.clone(),
        }
    }
    pub fn substitute(&self, range: &VarRange<F>, id: usize) -> Self {
        match self.substitute0(range, id) {
            Some(e) => e,
            None => self.clone(),
        }
    }
    fn substitute0(&self, range: &VarRange<F>, id: usize) -> Option<Self> {
        fn rc<T: Clone>(op: Option<T>, default: &Rc<T>) -> Rc<T> {
            match op {
                Some(e) => Rc::new(e),
                None => default.clone(),
            }
        }
        match self {
            Var(_) => Some(range.eval(id)),
            Add(e1, e2) => match (e1.substitute0(range, id), e2.substitute0(range, id)) {
                (None, None) => None,
                (f1, f2) => Some(Add(rc(f1, e1), rc(f2, e2))),
            },
            Mul(e1, e2) => match (e1.substitute0(range, id), e2.substitute0(range, id)) {
                (None, None) => None,
                (f1, f2) => Some(Mul(rc(f1, e1), rc(f2, e2))),
            },
            Neg(e) => e.substitute0(range, id).map(|x| Neg(Rc::new(x))),
            Unlab(e) => e.substitute0(range, id).map(|x| Unlab(Rc::new(x))),
            Named(e, name, latex) => e
                .substitute0(range, id)
                .map(|x| Named(Rc::new(x), name.clone(), *latex)),
            FromFunction(_, _)
            | FromIndicator(_, _)
            | Flag(_, _)
            | Unknown
            | Num(_)
            | Zero
            | One => None,
        }
    }
}

impl<N, F: Flag> Expr<N, F> {
    pub fn map<Fun, M>(&self, f: &Fun) -> Expr<M, F>
    where
        Fun: Fn(&N) -> M,
    {
        let rec = |e: &Self| Rc::new(e.map(f));

        match self {
            Add(e1, e2) => Add(rec(e1), rec(e2)),
            Mul(e1, e2) => Mul(rec(e1), rec(e2)),
            Neg(e) => Neg(rec(e)),
            Unlab(e) => Unlab(rec(e)),
            Named(e, name, latex) => Named(rec(e), name.clone(), *latex),
            FromFunction(_g, b) => FromFunction(Rc::new(|_, _| unimplemented!()), *b), // Fixme
            FromIndicator(g, b) => FromIndicator(*g, *b),
            Var(i) => Var(*i),
            Flag(id, b) => Flag(*id, *b),
            Unknown => Unknown,
            Num(n) => Num(Rc::new(f(n))),
            Zero => Zero,
            One => One,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::Graph;
    #[test]
    fn test_eval_expr() {
        type V = QFlag<i64, Graph>;
        let basis = Basis::new(4);
        let flag1: V = basis.flag_from_id(3);
        let flag2: V = basis.qflag_from_coeff(|g, _| g.edges().count() as i64);
        let flag3: V = basis.qflag_from_indicator(|g, _| g.connected());
        let result = flag1 + (flag2 * 3) - flag3;
        let result2 = result.expr.eval();
        assert_eq!(result, result2);

        let t = Type::new(2, 1);
        let b = Basis::new(3).with_type(t);
        let flag: V = b.flag_from_id(1);
        let res = ((flag.clone() * 3) * -flag).untype();
        let res2 = res.expr.eval();
        assert_eq!(res, res2)
    }
}
