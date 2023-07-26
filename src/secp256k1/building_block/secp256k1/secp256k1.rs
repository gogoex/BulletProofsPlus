use crate::secp256k1::building_block::{
  field::prime_field::PrimeField,
  secp256k1::{
    affine_point::AffinePoint,
    equation::Equation,
  },
};
use std::rc::Rc;
use num_bigint::BigUint;

#[derive(Debug, Clone)]
pub struct Secp256k1 {
  pub f: PrimeField,    // base prime field
  pub f_n: PrimeField,  // field of order n for convenience
  pub n: BigUint,  // order of g
  pub eq: Equation,
}

impl Secp256k1 {
  pub fn new() -> Self {
    // base prime field
    let base_field_order = BigUint::parse_bytes(b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap();
    let f = PrimeField::new(&base_field_order);

    // order of the base point
    let n = BigUint::parse_bytes(b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap();
    let f_n = PrimeField::new(&n);


    let a1 = f.elem(&0u8);
    let a2 = f.elem(&0u8);
    let a3 = f.elem(&0u8);
    let a4 = f.elem(&0u8);
    let a6 = f.elem(&7u8);
    let eq = Equation::new(a1, a2, a3, a4, a6);

    Secp256k1 {
      f,
      f_n,
      n,
      eq,
    }
  }

  pub fn g(&self) -> AffinePoint {
    // TODO compute only once at the first call
    let x = BigUint::parse_bytes(b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap();
    let y = BigUint::parse_bytes(b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap();
    let x = &self.f.elem(&x);
    let y = &self.f.elem(&y);
    let g = AffinePoint::new(&Rc::new(self.clone()), x, y);
    g
  }
}
