use crate::bls12_381::building_block::{
  scalar::prime_field_elem::PrimeFieldElem,
  point::point::Point,
  zero::Zero2,
};

pub struct MulVec {
  scalars: Vec<PrimeFieldElem>,
  points: Vec<Point>,
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
  ) -> Point {
    if self.scalars.len() != self.points.len() {
      panic!("mulvec: lengths of scalars and points must match");
    }
    let mut sum = Point::zero();

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

  pub fn add_point(&mut self, point: &Point) {
    self.points.push(point.clone());
  }

  pub fn add_points(&mut self, points: &Vec<Point>) {
    for point in points {
      self.points.push(point.clone());
    }
  }
}
