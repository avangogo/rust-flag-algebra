use approx::assert_relative_eq;
use flag_algebra::{
    tools::{SdpaCertificate, CERTIFICATE_FILE},
    *,
};
use flags::Graph;

fn make_goodman_problem() -> Problem<f64, Graph> {
    let basis = Basis::new(3);

    Problem {
        ineqs: vec![total_sum_is_one(basis), flags_are_nonnegative(basis)],
        cs: basis.all_cs(),
        obj: flag(&Graph::clique(3)) + flag(&Graph::empty(3)),
    }
}

fn make_turan_problem(forbidden: &Graph, n: usize) -> Problem<f64, Graph> {
    let basis = Basis::new(n);

    Problem {
        ineqs: vec![
            total_sum_is_one(basis),
            flags_are_nonnegative(basis),
            flag(forbidden).expand(basis).equal(0.0),
        ],
        cs: basis.all_cs(),
        obj: -flag(&Graph::clique(2)).expand(basis),
    }
}

#[test]
pub fn solve_and_read_certificate() {
    let problem = make_goodman_problem();

    let temp_dir = tempfile::tempdir().unwrap();
    let file = format!("{}/goodman_e2e", temp_dir.path().to_str().unwrap());

    problem.write_sdpa(&file).unwrap();
    let result = problem.run_csdp(&file, None, false).unwrap();
    assert_eq!(result, 0.25);

    let certificate = SdpaCertificate::load(CERTIFICATE_FILE).unwrap();

    assert_relative_eq!(
        *certificate.y.as_slice(),
        [0.125, 0.375, 0.375, 0.125],
        epsilon = 0.00001
    );
    assert_eq!(certificate.coeffs.len(), 18);
}

#[test]
pub fn turan_with_csdp() {
    let triangle = Graph::clique(3);
    let problem = make_turan_problem(&triangle, 3);

    let temp_dir = tempfile::tempdir().unwrap();
    let file = format!("{}/es_e2e", temp_dir.path().to_str().unwrap());
    assert_eq!(problem.solve_csdp(&file).unwrap(), -0.5);
}
