use crate::secp256k1::building_block::{
  field::prime_field_elem::PrimeFieldElem,
  secp256k1::affine_point::AffinePoint,
  secp256k1::secp256k1::Secp256k1,
  zero::Zero,
};
use std::rc::Rc;

pub struct MulVec {
  scalars: Vec<PrimeFieldElem>,
  points: Vec<AffinePoint>,
}

impl MulVec {
  pub fn new() -> Self {
    MulVec {
      scalars: vec![],
      points: vec![],
    }
  }

  pub fn calculate(
    &self,
    curve: &Rc<Secp256k1>,
  ) -> AffinePoint {
    if self.scalars.len() != self.points.len() {
      panic!("mulvec: lengths of scalars and points must match");
    }
    let mut sum = curve.g().zero();

    for i in 0..self.scalars.len() {
      let p = &self.points[i] * &self.scalars[i];
      sum = &sum + &p;
    }
    return sum;
  }

  pub fn add_scalar(&mut self, scalar: &PrimeFieldElem) {
    self.scalars.push(scalar.clone());
  }

  pub fn add_scalars(&mut self, scalars: &Vec<PrimeFieldElem>) {
    for scalar in scalars {
        self.scalars.push(scalar.clone());
    }
  }

  pub fn add_point(&mut self, point: &AffinePoint) {
    self.points.push(point.clone());
  }

  pub fn add_points(&mut self, points: &Vec<AffinePoint>) {
    for point in points {
      self.points.push(point.clone());
    }
  }
}
