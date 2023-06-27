extern crate flag_algebra;

use crate::flags::*;
use flag_algebra::*;

fn identity1<F>(t: Type<F>, n: usize)
where
    F: Flag,
{
    let b = Basis::new(n).with_type(t);
    let x: QFlag<i64, _> = b.random();
    let y = b.random();
    assert_eq!(&(&x - &y) * &(x), (&(&x * &x) - &(&y * &x)));
}

fn commute<F>(t: Type<F>, n1: usize, n2: usize)
where
    F: Flag,
{
    let b1 = Basis::new(n1).with_type(t);
    let b2 = Basis::new(n2).with_type(t);
    let x: QFlag<i64, _> = b1.random();
    let y: QFlag<i64, _> = b2.random();
    assert_eq!(&x * &y, &y * &x)
}

fn assoc<F>(t: Type<F>, n1: usize, n2: usize, n3: usize)
where
    F: Flag,
{
    let b1 = Basis::new(n1).with_type(t);
    let b2 = Basis::new(n2).with_type(t);
    let b3 = Basis::new(n3).with_type(t);
    let x: QFlag<i64, _> = b1.random();
    let y: QFlag<i64, _> = b2.random();
    let z: QFlag<i64, _> = b3.random();
    assert_eq!(&(&x * &y) * &z, &x * &(&y * &z))
}

fn transitive<F>(t: Type<F>, n1: usize, n2: usize, n3: usize)
where
    F: Flag,
{
    let b1 = Basis::new(n1).with_type(t);
    let b2 = Basis::new(n2).with_type(t);
    let b3 = Basis::new(n3).with_type(t);
    let x: QFlag<i64, _> = b1.random();
    assert_eq!(x.expand(b2).expand(b3), x.expand(b3))
}

#[test]
pub fn identities() {
    let t1 = Type::new(1, 0);
    identity1::<Graph>(Type::empty(), 2);
    identity1::<Graph>(t1, 2);
}

#[test]
pub fn mul_commutativity() {
    commute::<Graph>(Type::new(1, 0), 3, 2);
    commute::<Graph>(Type::empty(), 2, 3);
    commute::<OrientedGraph>(Type::new(1, 0), 3, 2);
}

#[test]
pub fn mul_associativity() {
    assoc::<Graph>(Type::new(1, 0), 2, 2, 2);
    assoc::<Graph>(Type::empty(), 2, 1, 2);
    assoc::<OrientedGraph>(Type::new(1, 0), 2, 1, 2);
}

#[test]
pub fn untype_transitivity() {
    transitive::<Graph>(Type::new(1, 0), 3, 4, 5);
    transitive::<Graph>(Type::empty(), 2, 3, 5);
    transitive::<OrientedGraph>(Type::new(1, 0), 3, 4, 5);
}
