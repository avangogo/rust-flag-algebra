//! WIP tools to minimize certificates

use crate::flag::Flag;
use crate::problem::sdp::*;
use crate::problem::sdpa;
use crate::tools::*;
use log::*;
use sprs::CsMat;

/// A higher level object to try several selectors on a problem
#[derive(Debug, Clone)]
pub struct FlagSolver<F: Flag> {
    pb: Problem<f64, F>,
    name: String,
    pub optimal_value: Option<f64>,
    /// If there is an optimal value, this correspond to a selector giving it
    select: Selector,
    /// The selector corresponding to the file name.spda
    select_sdpa_file: Option<Selector>,
    /// Selector corresponding to the file certificate
    select_certificate_file: Option<Selector>,
    protected: Vec<bool>,
}

impl<F> FlagSolver<F>
where
    F: Flag,
{
    /// Create a new instance that owns a problem
    pub fn new<S: Into<String>>(pb: Problem<f64, F>, name: S) -> Self {
        Self {
            protected: vec![false; pb.ineqs.len()],
            select: Selector::new(&pb), // Current best selector
            name: name.into(),
            pb,
            optimal_value: None,
            select_sdpa_file: None,
            select_certificate_file: None,
        }
    }
    /// Prevents inequality `i` to be ruled out
    pub fn protect(mut self, i: usize) -> Self {
        self.protected[i] = true;
        self
    }
    fn write_sdpa(&mut self, select: Selector) {
        self.pb
            .view(&select)
            .write_sdpa(&self.name)
            .expect("Cannot write sdpa file");
        self.select_sdpa_file = Some(select);
    }
    fn run_csdp(&mut self) -> Result<f64, ()> {
        match self.select_sdpa_file.take() {
            Some(select) => {
                self.select_certificate_file = Some(select);
                match self.pb.run_csdp(&self.name, None, false) {
                    Ok(v) => Ok(v),
                    Err(sdpa::Error::SdpNotSolved(_)) => Err(()),
                    Err(e) => panic!("Failed to run csdp {e}"),
                }
            }
            None => panic!("No problem to solve"),
        }
    }
    pub fn init(&mut self) {
        self.write_sdpa(self.select.clone());
        self.optimal_value = Some(self.run_csdp().expect("Cannot find initial solution"))
    }
    fn use_certificate(&mut self) {
        self.select = self
            .select
            .refine_with_certificate(&self.load_certificate(), &self.protected)
    }
    pub fn run(&mut self, select: Selector) -> Result<(), String> {
        self.write_sdpa(select.clone());
        if let Ok(v) = self.run_csdp() {
            if match self.optimal_value {
                None => {
                    self.optimal_value = Some(v);
                    true
                }
                Some(v0) => (v - v0).abs() < 1e-6,
            } {
                self.select = select;
                info!("New certificate of weight {:?}", self.select.weight());
                return Ok(());
            } else {
                trace!("Selector rejected");
                return Err("Selector rejected".into());
            }
        };
        Err("run_csdp failed".into())
    }
    fn load_certificate(&self) -> Certificate<f64> {
        if let Some(ref select) = self.select_certificate_file {
            Certificate::from_file_select(&self.pb.view(select), "certificate").unwrap()
        } else {
            panic!("No certificate yet")
        }
    }
    /// Try to minimize using the certificate
    pub fn minimize_certificate(&mut self) {
        self.write_sdpa(self.select.clone());
        let _ = self.run_csdp().expect("Cannot find the same value");
        let new_select = self
            .select
            .refine_with_certificate(&self.load_certificate(), &self.protected);
        self.run(new_select).expect("Certificate too simplified")
    }
    /// Try to remove each `cauchy_schwarz` to see if they are useful
    pub fn cs_elim(&mut self) {
        info!("Try to remove Cauchy-Schwarz");
        for &i in self.select.cs_vec().iter().rev() {
            if let Ok(select) = self.select.remove_cs(i) {
                let _ = self.run(select);
            } else {
                unreachable!()
            }
        }
    }
    pub fn thin_cs_elim(&mut self) {
        info!("Try to refine Cauchy-Schwarz");
        for &id in &self.select.cs_vec() {
            let dim = self.select.cs_dim(id);
            for i in 0..dim {
                let mut mat = CsMat::zero((dim, dim));
                mat.insert(i, i, 1.);
                if let Ok(select) = self.select.restrict_cs(id, mat) {
                    let _ = self.run(select);
                } else {
                    unreachable!()
                }
            }
        }
    }
    pub fn ineqs_elim(&mut self) {
        info!("Try to remove inequalities");
        for i in 0..self.select.ineqs.len() {
            if !self.protected[i] {
                let mut select = self.select.clone();
                select.ineqs[i].clear();
                let _ = self.run(select);
            }
        }
    }
    pub fn thin_ineqs_elim(&mut self) {
        info!("Try to refine inequalities");
        for i in (0..self.select.ineqs.len()).rev() {
            if !self.protected[i] {
                for j in 0..self.pb.ineqs[i].data.len() {
                    if let Ok(posj) = self.select.ineqs[i].binary_search(&j) {
                        let mut select = self.select.clone();
                        let _ = select.ineqs[i].remove(posj);
                        if let Ok(()) = self.run(select) {
                            self.use_certificate();
                        }
                    }
                }
            }
        }
    }
    pub fn print_report(&self)
    where
        F: Draw,
    {
        assert_eq!(Some(self.select.clone()), self.select_certificate_file);
        let cert = self.load_certificate();
        print_report(&self.pb.view(&self.select), &cert, "report").unwrap()
    }
    // Strategies

    pub fn minimize_once(&mut self)
    where
        F: Draw,
    {
        if self.optimal_value.is_none() {
            self.init();
        }
        self.print_report()
    }
    pub fn minimize(&mut self)
    where
        F: Draw,
    {
        if self.optimal_value.is_none() {
            self.init();
        }
        self.cs_elim();
        self.ineqs_elim();
        self.thin_ineqs_elim();
        self.write_sdpa(self.select.clone());
        let _ = self.run_csdp().unwrap();
        self.print_report()
    }
    pub fn minimize2(&mut self)
    where
        F: Draw,
    {
        if self.optimal_value.is_none() {
            self.init()
        };
        self.cs_elim();
        self.thin_cs_elim();
        self.ineqs_elim();
        self.minimize_certificate();
        self.write_sdpa(self.select.clone());
        let _ = self.run_csdp().unwrap();
        self.print_report()
    }
    pub fn minimize3(&mut self)
    where
        F: Draw,
    {
        info!("Method 3");
        if self.optimal_value.is_none() {
            self.init()
        };
        self.cs_elim();
        //self.thin_cs_elim();
        //self.minimize_certificate();
        self.ineqs_elim();
        self.thin_cs_elim();
        self.thin_ineqs_elim();
        self.write_sdpa(self.select.clone());
        let _ = self.run_csdp().unwrap();
        println!("{:?}", self.optimal_value);
        self.print_report()
    }
}
