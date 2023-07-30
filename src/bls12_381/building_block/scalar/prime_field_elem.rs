use mcl_rust::*;
use crate::bls12_381::building_block::zero::Zero;
use std::{
  fmt,
  ops::{Add, Sub, Mul, Div, Neg, Deref},
};
use bitvec::{
  prelude::Lsb0,
  view::BitView,
};

#[derive(Clone)]
pub struct PrimeFieldElem {
  pub e: Fr,
}

impl fmt::Debug for PrimeFieldElem {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{:?}", self.e.get_str(16))
  }
}

impl Zero<PrimeFieldElem> for PrimeFieldElem {
  fn zero(&self) -> Self {
    PrimeFieldElem {
      e: Fr::zero(),
    }
  }

  fn is_zero(&self) -> bool {
    self.e.is_zero()
  }
}

impl PartialEq for PrimeFieldElem {
  fn eq(&self, other: &Self) -> bool {
    self.e == other.e
  }
}

impl Eq for PrimeFieldElem {}

impl Deref for PrimeFieldElem {
  type Target = Fr;

  fn deref(&self) -> &Self::Target {
    &self.e
  }
}

macro_rules! impl_add {
  ($rhs: ty, $target: ty) => {
    impl<'a> Add<$rhs> for $target {
      type Output = PrimeFieldElem;

      fn add(self, rhs: $rhs) -> Self::Output {
        PrimeFieldElem {
          e: self.e.add(&rhs.e)
        }
      }
    }
  };
}
impl_add!(PrimeFieldElem, &PrimeFieldElem);
impl_add!(&PrimeFieldElem, &PrimeFieldElem);
impl_add!(&PrimeFieldElem, PrimeFieldElem);
impl_add!(PrimeFieldElem, PrimeFieldElem);

macro_rules! impl_sub {
  ($rhs: ty, $target: ty) => {
    impl<'a> Sub<$rhs> for $target {
      type Output = PrimeFieldElem;

      fn sub(self, rhs: $rhs) -> Self::Output {
        PrimeFieldElem {
          e: self.e.sub(&rhs.e)
        }
      }
    }
  };
}
impl_sub!(PrimeFieldElem, &PrimeFieldElem);
impl_sub!(&PrimeFieldElem, &PrimeFieldElem);
impl_sub!(&PrimeFieldElem, PrimeFieldElem);
impl_sub!(PrimeFieldElem, PrimeFieldElem);

macro_rules! impl_mul {
  ($rhs: ty, $target: ty) => {
    impl<'a> Mul<$rhs> for $target {
      type Output = PrimeFieldElem;

      fn mul(self, rhs: $rhs) -> Self::Output {
        PrimeFieldElem {
          e: self.e.mul(&rhs.e)
        }
      }
    }
  };
}
impl_mul!(PrimeFieldElem, &PrimeFieldElem);
impl_mul!(&PrimeFieldElem, &PrimeFieldElem);
impl_mul!(&PrimeFieldElem, PrimeFieldElem);
impl_mul!(PrimeFieldElem, PrimeFieldElem);

macro_rules! impl_div {
  ($rhs: ty, $target: ty) => {
    impl<'a> Div<$rhs> for $target {
      type Output = PrimeFieldElem;

      fn div(self, rhs: $rhs) -> Self::Output {
        PrimeFieldElem {
          e: self.e.div(&rhs.e)
        }
      }
    }
  };
}
impl_div!(PrimeFieldElem, &PrimeFieldElem);
impl_div!(&PrimeFieldElem, &PrimeFieldElem);
impl_div!(&PrimeFieldElem, PrimeFieldElem);
impl_div!(PrimeFieldElem, PrimeFieldElem);

// macro_rules! impl_bit_and {
//   ($rhs: ty, $target: ty) => {
//     impl<'a> BitAnd<$rhs> for $target {
//       type Output = PrimeFieldElem;

//       fn bitand(self, rhs: $rhs) -> Self::Output {
//         let res = &self.e & rhs.e.clone();
//         PrimeFieldElem {
//           f: self.f.clone(),
//           e: res,
//         }
//       }
//     }
//   }
// }
// impl_bit_and!(PrimeFieldElem, &PrimeFieldElem);
// impl_bit_and!(&PrimeFieldElem, &PrimeFieldElem);
// impl_bit_and!(&PrimeFieldElem, PrimeFieldElem);
// impl_bit_and!(PrimeFieldElem, PrimeFieldElem);

// macro_rules! impl_shr_assign {
//   ($rhs: ty, $target: ty) => {
//     impl ShrAssign<$rhs> for $target {
//       fn shr_assign(&mut self, rhs: $rhs) {
//         let n = rhs.to_u64().unwrap();
//         self.e >>= n;
//       }
//     }
//   }
// }
// impl_shr_assign!(BigUint, PrimeFieldElem);
// // impl_shr_assign!(BigUint, &PrimeFieldElem);
// impl_shr_assign!(&BigUint, PrimeFieldElem);
// // impl_shr_assign!(&BigUint, &PrimeFieldElem);

impl Neg for PrimeFieldElem {
  type Output = Self;

  fn neg(self) -> Self::Output {
    PrimeFieldElem {
      e: Fr::zero().sub(&self.e)
    }
  }
}

impl<'a> Neg for &'a PrimeFieldElem {
  type Output = PrimeFieldElem;

  fn neg(self) -> Self::Output {
    let mut e = Fr::zero();
    Fr::neg(&mut e, &self.e);

    PrimeFieldElem {
      e: Fr::zero().sub(&self.e)
    }
  }
}

impl From<Fr> for PrimeFieldElem {
  fn from(value: Fr) -> Self {
      PrimeFieldElem { e: value }
  }
}

impl PrimeFieldElem {
  pub fn init() {
    let b = init(CurveType::BLS12_381);
    if !b {
        panic!("Initializing mcl scalar failed");
    }
  }

  pub fn new(n: i32) -> Self {
    let mut e = Fr::zero();
    e.set_int(n);
    PrimeFieldElem { e }
  }

  pub fn one() -> Self {
    PrimeFieldElem::new(1)
  }

  // calculate w/ binary method
  pub fn pow(&self, rhs: u32) -> PrimeFieldElem {
    let rhs_le_bytes = rhs.to_le_bytes();

    let mut prod = PrimeFieldElem::one().e;

    let mut bit_value = self.e.clone();
    let rhs_in_bits = rhs_le_bytes.view_bits::<Lsb0>();

    for bit in rhs_in_bits {
      if bit == true {
        prod = &prod * &bit_value;
      }
      bit_value = &bit_value * &bit_value
    }

    PrimeFieldElem { e: prod }
  }

  pub fn sq(&self) -> PrimeFieldElem {
    let mut e = Fr::zero();
    Fr::sqr(&mut e, &self.e);

    PrimeFieldElem {
      e,
    }
  }

  pub fn inv(&self) -> PrimeFieldElem {
    let mut e = Fr::zero();
    Fr::inv(&mut e, &self.e);

    PrimeFieldElem {
      e,
    }
  }

  // 2. Compute 1/(e_1...e_k) and 1/e_k, ..., 1/e_1
  pub fn batch_invert(xs: &Vec<PrimeFieldElem>) -> (PrimeFieldElem, Vec<PrimeFieldElem>) {
    let mut prod = PrimeFieldElem::one().e;
    let mut inv_xs = vec![];

    for x in xs {
      inv_xs.push(x.inv());
      prod = prod.mul(&x.inv().e);
    }
    (prod.into(), inv_xs)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use std::sync::Once;

  static INIT: Once = Once::new();

  pub fn initialize() {
    INIT.call_once(|| {
      PrimeFieldElem::init();
    });
  }

  #[test]
  fn sub() {
    initialize();
    let a = PrimeFieldElem::new(9);
    let b = PrimeFieldElem::new(2);
    let c = a - b;

    let exp = PrimeFieldElem::new(7);
    assert_eq!(c, exp);
  }

  #[test]
  fn sub_eq_val() {
    initialize();
    let a = PrimeFieldElem::new(9);
    let b = PrimeFieldElem::new(9);
    let c = a - b;

    let exp = PrimeFieldElem::new(0);
    assert_eq!(c, exp);
  }

  #[test]
  fn mul() {
    initialize();
    let a = PrimeFieldElem::new(2);
    let b = PrimeFieldElem::new(5);
    let c = a * b;

    let exp = PrimeFieldElem::new(10);
    assert_eq!(c, exp);
  }

  #[test]
  fn div() {
    initialize();
    let a = PrimeFieldElem::new(10);
    let b = PrimeFieldElem::new(2);
    let c = a.div(&b);

    let exp = PrimeFieldElem::new(5);
    assert_eq!(c, exp);
  }

  #[test]
  fn inv() {
    initialize();
    let a = PrimeFieldElem::new(5);
    let inv_a = a.inv();
    let c = a * inv_a;

    let exp = PrimeFieldElem::new(1);
    assert_eq!(c, exp);
  }

  #[test]
  fn neg() {
    initialize();
    let a = PrimeFieldElem::new(5);
    let neg_a = a.clone().neg();
    let c = &a + neg_a;
    assert!(c.is_zero());
  }

  #[test]
  fn pow() {
    initialize();
    {
      let a = PrimeFieldElem::new(3);
      let c = a.pow(10);

      let exp = PrimeFieldElem::new(59049);
      assert_eq!(c, exp);
    }
    {
      let a = PrimeFieldElem::new(3);
      let c = a.pow(4);

      let exp = PrimeFieldElem::new(81);
      assert_eq!(c, exp);
    }
  }

  #[test]
  fn sq() {
    initialize();
    let a = PrimeFieldElem::new(9);
    let a_sq = a.sq();

    let exp = PrimeFieldElem::new(81);
    assert_eq!(a_sq, exp);
  }
}
