#![allow(non_snake_case)]
use crate::{
    publickey::PublicKey,
    bls12_381::building_block::{
        point::point::Point,
        scalar::prime_field_elem::PrimeFieldElem,
    },
};

/**
 * Range Prover which contains witness
 */
pub struct RangeProver {
    pub v_vec: Vec<u64>,
    pub gamma_vec: Vec<PrimeFieldElem>,
    pub commitment_vec: Vec<Point>,
}

impl RangeProver {
    pub fn new() -> Self {
        RangeProver {
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
        let commitment = pk.commitment(
            &PrimeFieldElem::new(v as i32),
            &gamma
        );
        println!("V[{}]={:?}", &self.v_vec.len() - 1, commitment);
        self.commitment_vec.push(commitment);
    }
}