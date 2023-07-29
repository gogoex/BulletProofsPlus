use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::Identity;

pub struct MulVec {
  scalars: Vec<Scalar>,
  points: Vec<RistrettoPoint>,
}

impl MulVec {
  pub fn new() -> Self {
    MulVec {
      scalars: vec![],
      points: vec![],
    }
  }

  fn mulvec(
    scalars: &Vec<Scalar>,
    points: &Vec<RistrettoPoint>,
  ) -> RistrettoPoint {
    if scalars.len() != points.len() {
      panic!("mulvec: lengths of scalars and points must match");
    }
    let mut sum = CompressedRistretto::identity().decompress().unwrap();

    for i in 0..scalars.len() {
      let p = &points[i] * &scalars[i];
      sum = &sum + &p;
    }
    return sum;
  }

  pub fn add_scalar(&mut self, scalar: &Scalar) {
    self.scalars.push(scalar.clone());
  }

  pub fn add_scalars(&mut self, scalars: &Vec<Scalar>) {
    for scalar in scalars {
        self.scalars.push(scalar.clone());
    }
  }

  pub fn add_point(&mut self, point: &RistrettoPoint) {
    self.points.push(point.clone());
  }

  pub fn add_points(&mut self, points: &Vec<RistrettoPoint>) {
    for point in points {
      self.points.push(point.clone());
    }
  }

  pub fn calculate(&self) -> RistrettoPoint {
    Self::mulvec(&self.scalars, &self.points)
  }
}
