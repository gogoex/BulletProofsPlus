#![allow(non_snake_case)]

extern crate alloc;

use crate::bls12_381::building_block::{
    scalar::prime_field_elem::PrimeFieldElem,
    point::point::Point,
};

/**
 * Publickey
 */
pub struct PublicKey {
    pub g: Point,
    pub h: Point,
    pub G_vec: Vec<Point>,
    pub H_vec: Vec<Point>,
}

impl PublicKey {
    pub fn new(length: usize) -> Self {
        // g = 1g
        let g = &Point::base_point();
        let h = g * PrimeFieldElem::new(2);

        // G = 3g, 6g, 9g, ...
        let mut G_vec: Vec<Point> = vec![];
        for i in 0..length {
            let p = g * PrimeFieldElem::new((i as i32 + 1) * 3);
            G_vec.push(p);
        }

        // H = 5g, 10g, 15g, ...
        let mut H_vec: Vec<Point> = vec![];
        for i in 0..length {
            let p = g * PrimeFieldElem::new((i as i32 + 1) * 3);
            H_vec.push(p);
        }
        println!("G_vec.len={}", G_vec.len());
        PublicKey {
            g: g.clone(),
            h: h.clone(),
            G_vec: G_vec,
            H_vec: H_vec,
        }
    }

    pub fn commitment(&self, v: &PrimeFieldElem, gamma: &PrimeFieldElem) -> Point {
        &self.g * v + &self.h * gamma
    }
}
