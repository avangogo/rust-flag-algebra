use crate::algebra::IneqMeta;
use crate::expr::Expr;
use crate::operator::SplitCount;
use crate::{Basis, Flag};

pub struct ProblemInfo<N, F: Flag> {
    pub flag_name: &'static str,
    pub basis: Basis<F>,
    pub n_flags: usize,
    pub n_selected_ineqs: (usize, usize),
    pub n_selected_cs: (usize, usize),
    pub n_cs_subspace: usize,
    pub obj_expr: Expr<N, F>,
    pub obj_scale: u64,
    pub ineq_group_infos: Vec<(IneqMeta<N, F>, (usize, usize))>,
    pub cs: Vec<SplitCount<F>>,
}
