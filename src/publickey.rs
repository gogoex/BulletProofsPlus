#![allow(non_snake_case)]

extern crate alloc;

use core::iter;
use crate::secp256k1::building_block::secp256k1::{
    affine_point::AffinePoint,
    secp256k1::Secp256k1,
};
use curve25519_dalek::traits::MultiscalarMul;

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
    //
    pub fn new(length: usize) -> Self {
        let curve = Secp256k1::new();

        let g = curve.g();
        let h = curve.g() * 2;

        let G_vec: Vec<AffinePoint> = curve.f.rand_elems(length, true);
        let H_vec: Vec<AffinePoint> = curve.f.rand_elems(length, true);

        PublicKey {
            g: g,
            h: h,
            G_vec: G_vec,
            H_vec: H_vec,
        }
    }

    pub fn commitment(&self, v: &Scalar, gamma: &Scalar) -> AffinePoint {
        self.g * v + self.h * gamma
    }

    pub fn vector_commitment(
        &self,
        a_vec: &Vec<Scalar>,
        b_vec: &Vec<Scalar>,
        out: &Scalar,
        gamma: &Scalar,
    ) -> RistrettoPoint {
        let scalars = a_vec
            .iter()
            .chain(b_vec.iter())
            .chain(iter::once(out))
            .chain(iter::once(gamma));
        let points = self.G_vec.iter()
            .chain(self.H_vec.iter())
            .chain(iter::once(&self.g))
            .chain(iter::once(&self.h));
            RistrettoPoint::multiscalar_mul(scalars, points)
    }
}