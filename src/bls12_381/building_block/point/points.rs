use crate::bls12_381::building_block::{
  scalar::{
    prime_field_elem::PrimeFieldElem,
    prime_field_elems::PrimeFieldElems,
  },
  point::point::Point,
  zero::Zero2,
};
use std::{
  fmt,
  ops::{Add, Mul, Deref},
};

pub struct Points {
  points: Vec<Point>,
}

impl Deref for Points {
  type Target = Vec<Point>;

  fn deref(&self) -> &Self::Target {
    &self.points
  }
}

impl Points {
  pub fn new(points: &Vec<Point>) -> Self {
    Points {
      points: points.clone(),
    }
  }

  pub fn sum(&self) -> Point {
    let mut sum = Point::zero();
    for p in &self.points {
      sum = sum + p;
    }
    sum
  }

  pub fn from(&self, idx: usize) -> Self {
    if idx >= self.len() {
      panic!("Index outside the range is specified");
    } else {
      let mut points = vec![];
      for i in idx..self.len() {
        points.push(self[i].clone());
      }
      Points::new(&points)
    }
  }

  pub fn to(&self, idx: usize) -> Self {
    if idx > self.points.len() {
      panic!("Index outside the range is specified");
    }
    let mut points = vec![];
    for i in 0..idx {
      points.push(self[i].clone());
    }
    Points::new(&points)
  }
}

impl fmt::Debug for Points {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      write!(f, "{{")?;
      for x in &self.points {
        write!(f, "{:?},", x)?;
      }
      write!(f, "}}")
  }
}

macro_rules! impl_add {
  ($rhs: ty, $target: ty) => {
    impl<'a> Add<$rhs> for $target {
      type Output = Points;

      fn add(self, rhs: $rhs) -> Self::Output {
        if self.len() != rhs.len() {
          panic!("Tried to add AffinePoints of diffrent length");
        }
        let mut points = vec![];
        for i in 0..self.len() {
          points.push(&self.points[i] + &rhs.points[i]);
        }
        Points::new(&points)
      }
    }
  };
}
impl_add!(Points, &Points);
impl_add!(&Points, &Points);
impl_add!(Points, Points);
impl_add!(&Points, Points);

macro_rules! impl_scalar_mul {
  ($rhs: ty, $target: ty) => {
    impl<'a> Mul<$rhs> for $target {
      type Output = Points;

      fn mul(self, rhs: $rhs) -> Self::Output {
        let mut points = vec![];
        for x in &self.points {
          points.push(x * rhs.clone())
        }
        Points::new(&points)
      }
    }
  };
}
impl_scalar_mul!(PrimeFieldElem, Points);
impl_scalar_mul!(&PrimeFieldElem, Points);
impl_scalar_mul!(PrimeFieldElem, &Points);
impl_scalar_mul!(&PrimeFieldElem, &Points);

macro_rules! impl_vec_mul {
  ($rhs: ty, $target: ty) => {
    impl<'a> Mul<$rhs> for $target {
      type Output = Points;

      fn mul(self, rhs: $rhs) -> Self::Output {
        if self.points.len() != rhs.len() {
          panic!("Tried to multiply PrimeFieldElems of different size to AffinePoints");
        }
        let mut points = vec![];
        for i in 0..self.points.len() {
          points.push(&self.points[i] * &rhs[i])
        }
        Points::new(&points)
      }
    }
  };
}
impl_vec_mul!(PrimeFieldElems, Points);
impl_vec_mul!(&PrimeFieldElems, Points);
impl_vec_mul!(PrimeFieldElems, &Points);
impl_vec_mul!(&PrimeFieldElems, &Points);

impl PartialEq for Points {
  fn eq(&self, rhs: &Self) -> bool {
    if self.points.len() != rhs.points.len() {
      false
    } else {
      for i in 0..self.points.len() {
        if self.points[i] != rhs.points[i] {
          return false;
        }
      }
      true
    }
  }
}

impl Eq for Points {}

#[cfg(test)]
mod tests {
  use crate::bls12_381::building_block::{
    point::{
      point::Point,
      points::Points,
    },
    arith::Arith,
  };

  #[test]
  fn test_from() {
    Arith::init();
    let g = &Point::base_point();
    let gs_vec = vec![
      g.clone(),
      g + g,
      g + g + g,
      g + g + g + g,
    ];
    let elems = Points::new(&gs_vec);
    {
      let res = &elems.from(0);
      assert_eq!(res.len(), 4);
      assert_eq!(res.as_slice(), gs_vec.as_slice());
    }
    {
      let res = &elems.from(1);
      assert_eq!(res.len(), 3);
      assert_eq!(&res[0], &gs_vec[1]);
      assert_eq!(&res[1], &gs_vec[2]);
      assert_eq!(&res[2], &gs_vec[3]);
    }
    {
      let res = &elems.from(2);
      assert_eq!(res.len(), 2);
      assert_eq!(&res[0], &gs_vec[2]);
      assert_eq!(&res[1], &gs_vec[3]);
    }
    {
      let res = &elems.from(3);
      assert_eq!(res.len(), 1);
      assert_eq!(&res[0], &gs_vec[3]);
    }
    // TODO test elem.from(4) and confirm it panics
  }

  #[test]
  fn test_to() {
    Arith::init();
    let g = &Point::base_point();
    let gs_vec = vec![
      g.clone(),
      g + g,
      g + g + g,
      g + g + g + g,
    ];
    let elems = Points::new(&gs_vec);
    {
      let res = &elems.to(0);
      assert_eq!(res.len(), 0);
    }
    {
      let res = &elems.to(1);
      assert_eq!(res.len(), 1);
      assert_eq!(&res[0], &gs_vec[0]);
    }
    {
      let res = &elems.to(2);
      assert_eq!(res.len(), 2);
      assert_eq!(&res[0], &gs_vec[0]);
      assert_eq!(&res[1], &gs_vec[1]);
    }
    {
      let res = &elems.to(3);
      assert_eq!(res.len(), 3);
      assert_eq!(&res[0], &gs_vec[0]);
      assert_eq!(&res[1], &gs_vec[1]);
      assert_eq!(&res[2], &gs_vec[2]);
    }
    {
      let res = &elems.to(4);
      assert_eq!(res.len(), 4);
      assert_eq!(&res[0], &gs_vec[0]);
      assert_eq!(&res[1], &gs_vec[1]);
      assert_eq!(&res[2], &gs_vec[2]);
      assert_eq!(&res[3], &gs_vec[3]);
    }
    // TODO test elem.to(5) and confirm it panics
  }

  #[test]
  fn test_sum() {
    Arith::init();
    let g = &Point::base_point();
    let gs_vec = vec![
      g.clone(),
      g + g,
    ];
    let elems = Points::new(&gs_vec);
    let act = &elems.sum();
    let exp = g + g + g;
    assert_eq!(act, &exp);
  }
}