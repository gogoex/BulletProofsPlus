#![allow(non_snake_case)]

extern crate alloc;

use crate::secp256k1::building_block::{
    field::prime_field_elem::PrimeFieldElem,
    secp256k1::{
        affine_point::AffinePoint,
        secp256k1::Secp256k1,
    },
};

/**
 * Publickey
 */
pub struct PublicKey {
    pub g: AffinePoint,
    pub h: AffinePoint,
    pub G_vec: Vec<AffinePoint>,
    pub H_vec: Vec<AffinePoint>,
}

impl PublicKey {
    pub fn new(length: usize) -> Self {
        let curve = Secp256k1::new();

        let g = curve.g();
        let h = curve.g() * curve.f.elem(&2u8);

        let mut G_vec: Vec<AffinePoint> = vec![];
        for _i in 0..length {
            let p = curve.g() * curve.f.rand_elem(true);
            G_vec.push(p);
        }
        let mut H_vec: Vec<AffinePoint> = vec![];
        for _i in 0..length {
            let p = curve.g() * curve.f.rand_elem(true);
            H_vec.push(p);
        }

        PublicKey {
            g: g,
            h: h,
            G_vec: G_vec,
            H_vec: H_vec,
        }
    }

    pub fn commitment(&self, v: &PrimeFieldElem, gamma: &PrimeFieldElem) -> AffinePoint {
        &self.g * v + &self.h * gamma
    }
}