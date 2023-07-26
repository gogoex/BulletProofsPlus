#![allow(non_snake_case)]
use crate::publickey::PublicKey;
use crate::secp256k1::building_block::{
    field::prime_field_elem::PrimeFieldElem,
    secp256k1::{
        affine_point::AffinePoint,
        secp256k1::Secp256k1,
    },
};
use std::rc::Rc;

/**
 * Range Prover which contains witness
 */
pub struct RangeProver {
    pub curve: Rc<Secp256k1>,
    pub(crate) v_vec: Vec<u64>,
    pub(crate) gamma_vec: Vec<PrimeFieldElem>,
    pub commitment_vec: Vec<AffinePoint>,
}

impl RangeProver {
    pub fn new(curve: Rc<Secp256k1>) -> Self {
        RangeProver {
            curve,
            v_vec: Vec::new(),
            gamma_vec: Vec::new(),
            commitment_vec: Vec::new(),
        }
    }

    pub fn commit(
        &mut self,
        pk: &PublicKey,
        v: u64,
        gamma: PrimeFieldElem,
    ) {
        self.v_vec.push(v);
        self.gamma_vec.push(gamma.clone());
        self.commitment_vec.push(
            pk.commitment(
                &PrimeFieldElem::new(&Rc::new(self.curve.f.clone()), &v),
                &gamma
            )
        );
    }
}