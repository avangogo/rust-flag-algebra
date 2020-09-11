use crate::flag::SubClass;
use crate::flags::{CGraph, Colored, Graph};
use svg::node::element::{Circle, Line, SVG};
use svg::node::Node;

/// Trait to draw flags in svg
pub trait Draw: Sized {
    fn draw_with_parameters<C>(&self, color: C, type_size: usize) -> SVG
    where
        C: FnMut(usize) -> usize;
    fn draw(&self) -> SVG {
        self.draw_with_parameters(|_| 0, 0)
    }
    fn draw_typed(&self, type_size: usize) -> SVG {
        self.draw_with_parameters(|_| 0, type_size)
    }
    fn to_svg_file(&self, filename: &str) {
        svg::save(filename, &self.draw()).unwrap()
    }
}

/// Inheritance
impl<F, A> Draw for SubClass<F, A>
where
    F: Draw,
{
    fn draw_with_parameters<C>(&self, color: C, type_size: usize) -> SVG
    where
        C: FnMut(usize) -> usize,
    {
        self.content.draw_with_parameters(color, type_size)
    }
}

impl<F, A> Draw for Colored<F, A>
where
    F: Draw,
{
    fn draw_with_parameters<C>(&self, mut color: C, type_size: usize) -> SVG
    where
        C: FnMut(usize) -> usize,
    {
        let n = self.color.len();
        assert!((0..n).all(|v| { color(v) == 0 }));
        //
        self.content
            .draw_with_parameters(|v| self.color[v] as usize, type_size)
    }
}

/// Particular implementations

fn coordinates(i: usize, n: usize) -> (f64, f64) {
    assert!(i < n);
    if n == 1 {
        (50., 50.)
    } else {
        let angle = 2. * std::f64::consts::PI * i as f64 / n as f64;
        (50. + 41. * angle.cos(), 50. + 41. * angle.sin())
    }
}

fn vertex(i: usize, n: usize) -> Circle {
    let (cx, cy) = coordinates(i, n);
    Circle::new()
        .set("r", 6.)
        .set("stroke", "black")
        .set("cx", cx)
        .set("cy", cy)
}

fn type_marker(i: usize, n: usize) -> Circle {
    let (cx, cy) = coordinates(i, n);
    Circle::new()
        .set("r", 8.5)
        .set("stroke", "black")
        .set("stroke-width", 1.5)
        .set("fill", "none")
        .set("cx", cx)
        .set("cy", cy)
}

fn line(i: usize, j: usize, n: usize) -> Line {
    let (x1, y1) = coordinates(i, n);
    let (x2, y2) = coordinates(j, n);
    Line::new()
        .set("x1", x1)
        .set("x2", x2)
        .set("y1", y1)
        .set("y2", y2)
        .set("stroke-width", 4)
}

fn color(c: usize) -> &'static str {
    ["black", "red", "blue"][c]
}

fn frame() -> SVG {
    SVG::new().set("viewBox", "0 0 100 100").set("width", 100)
}

fn add_vertex<N: Node>(node: &mut N, i: usize, n: usize, col: &'static str, in_type: bool) {
    node.append(vertex(i, n).set("fill", col));
    if in_type {
        node.append(type_marker(i, n));
    }
}

impl<E> Draw for CGraph<E> {
    fn draw_with_parameters<C>(&self, mut col: C, type_size: usize) -> SVG
    where
        C: FnMut(usize) -> usize,
    {
        let mut res = frame();
        let n = self.size;
        for u in 0..n {
            for v in 0..u {
                match self.edge(u, v) {
                    0 => (),
                    c => res = res.add(line(u, v, n).set("stroke", color(c as usize - 1))),
                }
            }
        }
        for v in 0..n {
            add_vertex(&mut res, v, n, color(col(v)), v < type_size)
        }
        res
    }
}

impl Draw for Graph {
    fn draw_with_parameters<C>(&self, mut col: C, type_size: usize) -> SVG
    where
        C: FnMut(usize) -> usize,
    {
        let mut res = frame();
        let n = self.size();
        for u in 0..n {
            for v in 0..u {
                if self.edge(u, v) {
                    res = res.add(line(u, v, n).set("stroke", "black"))
                }
            }
        }
        for v in 0..n {
            add_vertex(&mut res, v, n, color(col(v)), v < type_size)
        }
        res
    }
}
