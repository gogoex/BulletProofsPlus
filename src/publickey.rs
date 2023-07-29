#![allow(non_snake_case)]

extern crate alloc;

use crate::secp256k1::building_block::{
    field::prime_field_elem::PrimeFieldElem,
    secp256k1::{
        affine_point::AffinePoint,
        secp256k1::Secp256k1,
    },
};
use std::rc::Rc;

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
    pub fn new(curve: &Rc<Secp256k1>, length: usize) -> Self {
        // g = 1g
        let g = curve.g();
        let h = curve.g() * curve.f_n.elem(&2u8);

        // G = 3g, 6g, 9g, ...
        let mut G_vec: Vec<AffinePoint> = vec![];
        for i in 0..length {
            let p = curve.g() * curve.f_n.elem(&((i + 1) * 3));
            G_vec.push(p);
        }

        // H = 5g, 10g, 15g, ...
        let mut H_vec: Vec<AffinePoint> = vec![];
        for i in 0..length {
            let p = curve.g() * curve.f_n.elem(&((i + 1) * 3));
            H_vec.push(p);
        }
        println!("G_vec.len={}", G_vec.len());
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