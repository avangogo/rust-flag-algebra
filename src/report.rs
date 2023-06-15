use self::ndarray_linalg::{Eigh, UPLO};
use crate::algebra::*;
use crate::draw::Draw;
use crate::expr::Names;
use crate::flag::Flag;
use crate::operator::{Basis, Savable};
use crate::sdp::*;
use ndarray::ScalarOperand;
use num::{FromPrimitive, Num, Zero};

use sprs::CsMat;
use std::fmt::Display;
use std::fs::File;
use std::io::*;
use std::ops::AddAssign;
use std::path::PathBuf;
extern crate ndarray_linalg;
use ndarray::Array1;
use ndarray::Array2;

pub trait Html {
    fn print_html<W: Write>(&self, w: &mut W) -> Result<()>;

    const LATEX: bool = false;
    /// True if Mathjax need to be loaded

    fn html(&self, name: &str) -> Result<()> {
        let mut filename = PathBuf::from(name);
        let _ = filename.set_extension("html");
        let mut file = BufWriter::new(File::create(&filename)?);
        writeln!(file, "<!DOCTYPE html><html><head><title>{}</title>", name)?;
        if Self::LATEX {
            writeln!(file, "{}", MATHJAX)?;
        }
        writeln!(file, "<style>{}</style></head><body>", CSS)?;
        self.print_html(&mut file)?;
        writeln!(file, "</body>")
    }
}

impl<F: Draw> Html for F {
    fn print_html<W: Write>(&self, w: &mut W) -> Result<()> {
        writeln!(w, "{}", self.draw())
    }
}

const CSS: &str = "
:root {
    --color1: #a3bbdc;
    --color2: #dae4f1;
    --color3: #edf2f8;
    --color4: #ffe8d6;
    --darkcolor4: #ff971d;
}
h2 {
    color: var(--darkcolor4);
}
body {
    background-color: #133454;
}
.flags {
    background-color: var(--color3);
    margin: 5px 10px;
}
details[open] > summary {
    background-color: var(--color1);
}
details {
    background-color: var(--color2);
}
details[open] {
    padding-bottom: 5px;
    margin-bottom: 5px;
}
div.qflag_item {
    display: inline-block;
    margin: 10px;
    text-align: center;
}
.qflag_item > span {
    display: block;
}
div.inequality {
    border-left-style:solid;
    border-color: var(--darkcolor4);
    border-width: thick;
    background-color: var(--color4);
    margin: 10px 0px 10px;
    padding: 10px;
}
div.obj {
    border-radius: 25px;
    background-color: var(--color4);
    margin: 10px 0px 10px;
    padding: 10px;

}
svg.inline-flag {
    width: 60px;
    top: 20px;
    position: relative;
}
";

const MATHJAX: &str = "<script>
MathJax = {
  tex: {
    inlineMath: ['\\\\[', '\\\\]', ['$','$']],
  }
};
</script>
     <script type=\"text/javascript\" id=\"MathJax-script\" async
     src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\">
     </script>";

fn print_tab<W, F, G, L>(w: &mut W, tab: &[G], mut f: F, type_size: usize) -> Result<()>
where
    W: Write,
    F: FnMut(usize, &G) -> Option<L>,
    G: Draw,
    L: Display,
{
    writeln!(w, "<div class=\"flags\">")?;
    for (i, x) in tab.iter().enumerate() {
        if let Some(label) = f(i, x) {
            writeln!(w, "<div class=\"qflag_item\">")?;
            writeln!(
                w,
                "<span>{}</span>",
                x.draw_with_parameters(|_| 0, type_size)
            )?;
            writeln!(w, "<span>{}</span></div>", label)?;
        }
    }
    writeln!(w, "</div>")
}

impl<G: Draw, F, L> Html for (&[G], F)
where
    F: FnMut(&G) -> L + Clone,
    L: Display,
{
    const LATEX: bool = G::LATEX;

    fn print_html<W: Write>(&self, w: &mut W) -> Result<()> {
        let mut f = self.1.clone();
        print_tab(w, self.0, |_, x| Some(f(x)), 0)?;
        Ok(())
    }
}

impl<F> Html for Basis<F>
where
    F: Flag + Draw,
{
    const LATEX: bool = F::LATEX;

    fn print_html<W: Write>(&self, w: &mut W) -> Result<()> {
        print_tab(w, &self.get(), |i, _| Some(i), self.t.size)
    }
}

impl<F, N> Html for QFlag<N, F>
where
    F: Flag + Draw,
    N: Num + Clone + Display + FromPrimitive,
{
    const LATEX: bool = F::LATEX;

    fn print_html<W: Write>(&self, w: &mut W) -> Result<()> {
        let scale: N = N::from_u64(self.scale).unwrap();
        print_tab(
            w,
            &self.basis.get(),
            |i, _| {
                let val = self.data[i].clone();
                if val.is_zero() {
                    None
                } else {
                    Some(val / scale.clone())
                }
            },
            self.basis.t.size,
        )
    }
}

fn inlined_list<'a, I>(iter: I) -> String
where
    I: Iterator + Clone + 'a,
    I::Item: AsRef<str>,
{
    let mut res = String::new();
    let last = iter.clone().count() - 1;
    for (i, word) in iter.enumerate() {
        if i > 0 {
            if i == last {
                res += " and ";
            } else {
                res += ", ";
            }
        };
        res += word.as_ref();
    }
    res
}

impl<N, F> Html for Names<N, F>
where
    F: Flag + Draw,
    N: Num + Clone + Display + FromPrimitive,
{
    const LATEX: bool = F::LATEX;

    fn print_html<W: Write>(&self, w: &mut W) -> Result<()> {
        // Inlined defintions
        let mut inline_names = Vec::new();
        for (t, name) in self.types.iter() {
            let svg = Basis::<F>::new(t.size).get()[t.id].draw_typed(t.size);
            inline_names.push((name, svg))
        }
        for ((i, basis), name) in self.flags.iter() {
            let svg = basis.get()[*i].draw_typed(basis.t.size);
            inline_names.push((name, svg))
        }
        if !inline_names.is_empty() {
            let defs: Vec<_> = inline_names
                .into_iter()
                .map(|(name, svg)| format!("${}=$ {}", name, svg.set("class", "inline-flag")))
                .collect();
            writeln!(w, "<p>where {}.</p>", inlined_list(defs.iter()))?
        }
        // Long definitions
        let other_names: Vec<_> = self
            .sets
            .iter()
            .map(|(name, _, _)| name)
            .chain(self.functions.iter().map(|(name, _)| name))
            .map(|s| format!("${}$", s))
            .collect();
        if !other_names.is_empty() {
            writeln!(
                w,
                "<details><summary>See the definition of {}.</summary>",
                inlined_list(other_names.iter())
            )?;
            // Detailed content
            for (name, basis, flags) in &self.sets {
                writeln!(w, "<p>${}$ contains the following flags:</p>", name)?;
                writeln!(w, "<div class=\"flags\">")?;
                for flag in flags.iter() {
                    writeln!(w, "{}", flag.draw_typed(basis.t.size))?;
                }
                writeln!(w, "</div>")?;
            }
            for (name, qflag) in &self.functions {
                writeln!(w, "<p>${}$ takes the following values:</p>", name)?;
                qflag.print_html(w)?;
            }

            writeln!(w, "</details>")?;
        }
        Ok(())
    }
}

pub trait Approx: Clone + Zero {
    fn is_negligible(&self) -> bool;
    fn round(&self) -> Self {
        if self.is_negligible() {
            Self::zero()
        } else {
            self.clone()
        }
    }
}

impl Approx for f64 {
    fn is_negligible(&self) -> bool {
        self.abs() < 1e-8
    }
}

fn round<N, F>(vec: &QFlag<N, F>) -> QFlag<N, F>
where
    N: Approx + Zero + Clone,
{
    vec.map(N::round)
}

// <script type=text/x-mathjax-config> MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$']]}}); </script>

fn header<W: Write>(w: &mut W, title: &str) -> Result<()> {
    writeln!(
        w,
        "<!DOCTYPE html><html><head>
<meta charset=\"UTF-8\">
<title>{}</title>
{}
<style>
{}
</style></head><body>",
        title, MATHJAX, CSS,
    )
}

fn footer<W: Write>(w: &mut W) -> Result<()> {
    writeln!(w, "</body></html>")
}

pub fn print_report<N, F>(
    // /!\ select must be applied to pb only
    pb: &ProblemView<N, F>,
    cert: &Certificate<f64>,
    filename: &str,
) -> Result<()>
where
    F: Flag + Draw,
    N: Display + Num + Clone + FromPrimitive + AddAssign + ScalarOperand + From<f64> + Approx,
{
    let mut filename = PathBuf::from(filename);
    let _ = filename.set_extension("html");
    let mut w = BufWriter::new(File::create(&filename)?);
    // Header
    header(&mut w, "Report")?;
    // Objective
    writeln!(w, "<h2>Objective:</h2>")?;
    writeln!(w, "<div class=\"obj\">")?;
    writeln!(w, "<p>Minimize: {}</p>", pb.obj.expr)?;
    {
        let mut names = Names::new();
        writeln!(w, "<p>\\[{}\\]</p>", pb.obj.expr.latex(&mut names))?;
        if !names.is_empty() {
            names.print_html(&mut w)?;
        }
    }
    writeln!(
        w,
        "<details><summary>Flag expression of the objective.</summary><div>"
    )?;
    round(pb.obj).print_html(&mut w)?;
    writeln!(w, "</div></details>")?;
    writeln!(w, "</div>")?;

    writeln!(w, "<h2>Inequalities</h2>")?;

    for (block, ineqs) in pb.ineqs.iter().enumerate() {
        assert!(ineqs.len() > 0);
        writeln!(w, "<div class=\"inequality\">")?;
        writeln!(w, "<p>{:.6}</p>", ineqs.meta())?;
        {
            let mut names = Names::new();
            writeln!(w, "<p>\\[{}\\]</p>", ineqs.meta().latex(&mut names))?;
            names.print_html(&mut w)?;
        }
        write_diag_csmatrix(&mut w, &cert.x[block])?;
        writeln!(w, "<details><summary>Contribution.</summary>")?;
        let coeff: Vec<N> = cert
            .diag_coeffs(block, ineqs.len())
            .into_iter()
            .map(|x| x.into())
            .collect();
        let (lhs, bound) = condense(ineqs, &coeff);
        round(&lhs).print_html(&mut w)?;
        writeln!(w, "\n is at least {}.</details>", bound)?;
        writeln!(w, "</div>")?;
    }

    writeln!(w, "<h2>Cauchy-Schwarz</h2>")?;

    for (block, cs) in pb.cs.iter().enumerate() {
        writeln!(w, "<div class=\"inequality\">")?;
        writeln!(w, "<h4>{}</h4><p>", cs)?;
        svg::write(
            &mut w,
            &cs.1.unlabeling.basis.get()[cs.1.unlabeling.flag].draw(),
        )?;
        writeln!(w, "</p>")?;
        let mat = &cert.x[pb.ineqs.len() + block];
        write_csmatrix(&mut w, mat)?;
        //
        writeln!(w, "<details><summary>Eigenvalue decomposition</summary>")?;
        let a: Array2<f64> = mat.to_dense();
        let (eigenvalues, eigenvectors) = a.eigh(UPLO::Lower).unwrap();
        write_array1(&mut w, &eigenvalues)?;
        write_array2(&mut w, &eigenvectors)?;
        writeln!(w, "</details>")?;
        //
        writeln!(w, "<details><summary>Eigenvectors</summary>")?;
        for (lambda, vect) in &eigenvectors_qflags(&cs, &a) {
            if !lambda.is_negligible() {
                writeln!(w, "<h>λ = {}</h><p>", lambda)?;
                round(vect).print_html(&mut w)?;
                writeln!(w, "</p>")?;
            }
        }
        writeln!(w, "</details>")?;
        //
        writeln!(w, "<details><summary>Contribution.</summary>")?;
        round(&condense_cs(&cs, &a)).print_html(&mut w)?;
        writeln!(w, "</details>")?;
        writeln!(w, "</div>")?;
    }
    footer(&mut w)
}

fn write_csmatrix<W: Write, N: Display>(w: &mut W, mat: &CsMat<N>) -> Result<()> {
    writeln!(w, "<math><mfenced><mtable>")?;
    for i in 0..mat.rows() {
        writeln!(w, "<mtr>")?;
        for j in 0..mat.cols() {
            if let Some(x) = mat.get(i, j) {
                writeln!(w, "<mtd>{}</mtd>", x)?;
            } else {
                writeln!(w, "<mtd></mtd>")?;
            }
        }
        writeln!(w, "</mtr>")?;
    }
    writeln!(w, "</mtable></mfenced></math>")
}

fn write_diag_csmatrix<W: Write, N: Display>(w: &mut W, mat: &CsMat<N>) -> Result<()> {
    writeln!(w, "<math><mfenced><mtable><mtr>")?;
    for i in 0..mat.rows() {
        if let Some(x) = mat.get(i, i) {
            writeln!(w, "<mtd>{}</mtd>", x)?;
        } else {
            writeln!(w, "<mtd>0</mtd>")?;
        }
    }
    writeln!(w, "</mtr></mtable></mfenced></math>")
}

fn write_array2<W: Write, N: Display>(w: &mut W, mat: &Array2<N>) -> Result<()> {
    writeln!(w, "<math><mfenced><mtable>")?;
    for i in 0..mat.nrows() {
        writeln!(w, "<mtr>")?;
        for j in 0..mat.ncols() {
            writeln!(w, "<mtd>{}</mtd>", mat[(i, j)])?;
        }
        writeln!(w, "</mtr>")?;
    }
    writeln!(w, "</mtable></mfenced></math>")
}

fn write_array1<W: Write, N: Display>(w: &mut W, mat: &Array1<N>) -> Result<()> {
    writeln!(w, "<math><mfenced><mtable>")?;
    writeln!(w, "<mtr>")?;
    for x in mat {
        writeln!(w, "<mtd>{}</mtd>", x)?;
    }
    writeln!(w, "</mtr>")?;
    writeln!(w, "</mtable></mfenced></math>")
}

// For including quantum flags in latex reports
pub fn latexify<F, N>(qflag: &QFlag<N, F>, folder: &str)
where
    F: Flag + Draw,
    N: Num + Clone + Display + FromPrimitive,
{
    let mut path = PathBuf::from(folder);
    let scale: N = N::from_u64(qflag.scale).unwrap();
    for (i, (val, flag)) in qflag
        .data
        .iter()
        .zip(qflag.basis.get().into_iter())
        .enumerate()
    {
        if !val.is_zero() {
            let b = qflag.basis;
            let filename = format!("{}in{}t{}id{}", i, b.size, b.t.size, b.t.id);
            path.push(&filename);
            path.set_extension("svg");
            svg::save(&path, &flag.draw_typed(qflag.basis.t.size)).unwrap();
            let x = val.clone() / scale.clone();
            if x.is_one() {
                print!(" + \\flag{{{}}}", filename)
            } else {
                print!(" + {}\\cdot\\flag{{{}}}", val, filename)
            };
            assert!(path.pop());
        }
    }
}
