use mcl_rust::*;
use crate::bls12_381::building_block::{
  scalar::prime_field_elem::PrimeFieldElem,
  zero::Zero2,
};
use std::{
  fmt,
  ops::{Add, Sub, Mul, Neg},
};

#[derive(Clone)]
pub struct Point(G1);

impl Point {
  pub fn base_point() -> Self {
    let p = G1::from_str("1 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap();
    Point(p)
  }

  pub fn double(&self) -> Self {
    let mut r = G1::zero();
    G1::dbl(&mut r, &self.0);
    Point(r)
  }
}

impl fmt::Debug for Point {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{{ x: {:?}, y: {:?}, z: {:?} }}", &self.0.x, &self.0.y, self.0.z)
  }
}

macro_rules! impl_add {
  ($rhs: ty, $target: ty) => {
    impl Add<$rhs> for $target {
      type Output = Point;
      fn add(self, rhs: $rhs) -> Self::Output {
        let mut r = G1::zero();
        G1::add(&mut r, &self.0, &rhs.0);
        Point(r)
      }
    }
  }
}
impl_add!(Point, Point);
impl_add!(&Point, Point);
impl_add!(Point, &Point);
impl_add!(&Point, &Point);

macro_rules! impl_sub {
  ($rhs: ty, $target: ty) => {
    impl Sub<$rhs> for $target {
      type Output = Point;

      fn sub(self, rhs: $rhs) -> Self::Output {
        let mut r = G1::zero();
        G1::sub(&mut r, &self.0, &rhs.0);
        Point(r)
      }
    }
  }
}
impl_sub!(Point, Point);
impl_sub!(&Point, Point);
impl_sub!(Point, &Point);
impl_sub!(&Point, &Point);

macro_rules! impl_mul {
  ($rhs: ty, $target: ty) => {
    impl Mul<$rhs> for $target {
      type Output = Point;

      fn mul(self, rhs: $rhs) -> Self::Output {
        let mut r = G1::zero();
        G1::mul(&mut r, &self.0, &rhs.e);
        Point(r)
      }
    }
  }
}
impl_mul!(PrimeFieldElem, Point);
impl_mul!(&PrimeFieldElem, Point);
impl_mul!(PrimeFieldElem, &Point);
impl_mul!(&PrimeFieldElem, &Point);

macro_rules! impl_neg {
  ($target: ty) => {
    impl Neg for $target {
      type Output = Point;

      fn neg(self) -> Self::Output {
        let mut r = G1::zero();
        G1::neg(&mut r, &self.0);
        Point(r)
      }
    }
  }
}
impl_neg!(Point);
impl_neg!(&Point);

impl Zero2<Point> for Point {
  fn zero() -> Self {
    Point(G1::zero())
  }

  fn is_zero(&self) -> bool {
    self.0.is_zero()
  }
}

impl PartialEq for Point {
  fn eq(&self, other: &Self) -> bool {
    self.0 == other.0
  }
}

impl Eq for Point {}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::bls12_381::building_block::arith::Arith;

  #[test]
  fn scalar_mul() {
    Arith::init();
    let g = &Point::base_point();

    {
      let n = PrimeFieldElem::new(1);
      let act = g * &n;
      assert_eq!(&act, g);
    }
    {
      let n = PrimeFieldElem::new(2);
      let act = g * n;
      let exp = g + g;
      assert_eq!(act, exp);
    }
    {
      let n = PrimeFieldElem::new(3);
      let act = g * n;
      let exp = g + g + g;
      assert_eq!(act, exp);
    }
  }

  #[test]
  fn double() {
    Arith::init();
    let g = &Point::base_point();
    let dbl_g = g.double();
    let g2 = g + g;
    assert_eq!(dbl_g, g2);
  }

  #[test]
  fn neg() {
    Arith::init();
    let g = &Point::base_point();
    let neg_g = -g;
    let act = g + neg_g;
    assert_eq!(act, Point::zero());
  }

  #[test]
  fn sub() {
    Arith::init();
    let g = &Point::base_point();
    let act = g - g;
    assert_eq!(act, Point::zero());
  }

  #[test]
  fn add() {
    Arith::init();
    let g = &Point::base_point();
    let g2 = g + g;

    let act = g + g2;
    let exp = g + g + g;
    assert_eq!(act, exp);
  }
}