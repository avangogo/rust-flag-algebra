extern crate flag_algebra;

use flag_algebra::*;
use crate::flags::*;
use crate::operator::*;

fn identity1<F>(t: Type, n: usize)
    where F: Flag
{
    let b = Basis::<F>::make(n, t);
    let x : QFlag<i64,_> = b.random();
    let y = b.random();
    assert_eq!(&(&x - &y) * &(x),
               (&(&x * &x) - &((&y * &x))));
}

fn commute<F>(t: Type, n1: usize, n2: usize)
    where F: Flag
{
    let b1 = Basis::<F>::make(n1, t);
    let b2 = Basis::<F>::make(n2, t);
    let x : QFlag<i64,_> = b1.random();
    let y : QFlag<i64,_> = b2.random();
    assert_eq!(&x * &y, &y * &x)
}

fn assoc<F>(t: Type, n1: usize, n2: usize, n3: usize)
    where F: Flag
{
    let b1 = Basis::<F>::make(n1, t);
    let b2 = Basis::<F>::make(n2, t);
    let b3 = Basis::<F>::make(n3, t);
    let x : QFlag<i64,_> = b1.random();
    let y : QFlag<i64,_> = b2.random();
    let z : QFlag<i64,_> = b3.random();
    assert_eq!(&(&x * &y) * &z, &x * &(&y  * &z))
}

fn transitive<F>(t: Type, n1: usize, n2: usize, n3: usize)
    where F: Flag
{
    let b1 = Basis::<F>::make(n1, t);
    let b2 = Basis::<F>::make(n2, t);
    let b3 = Basis::<F>::make(n3, t);
    let x : QFlag<i64,_> = b1.random();
    assert_eq!( x.expand(b2).expand(b3), x.expand(b3) )
}

#[test]
pub fn identities() {
    let t1 = Type::new(1,0);
    identity1::<Graph>(Type::empty(), 2);
    identity1::<Graph>(t1, 2);
}

#[test]
pub fn mul_commutativity() {
    let t1 = Type::new(1,0);
    commute::<Graph>(t1, 3, 2);
    commute::<Graph>(Type::empty(), 2, 3);
    commute::<Digraph>(t1, 3, 2);
}

#[test]
pub fn mul_associativity() {
    let t1 = Type::new(1,0);
    assoc::<Graph>(t1, 2, 2, 2);
    assoc::<Graph>(Type::empty(), 2, 1, 2);
    assoc::<Digraph>(t1, 2, 1, 2);
}


#[test]
pub fn untype_transitivity() {
    let t1 = Type::new(1,0);
    transitive::<Graph>(t1, 3, 4, 5);
    transitive::<Graph>(Type::empty(), 2, 3, 5);
    transitive::<Digraph>(t1, 3, 4, 5);
}
