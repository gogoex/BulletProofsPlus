/**
 * This code is mostly copied from
 * https://github.com/dalek-cryptography/bulletproofs
 */

use crate::secp256k1::building_block::{
    field::prime_field_elem::PrimeFieldElem,
    secp256k1::secp256k1::Secp256k1,
};

use std::rc::Rc;

pub struct ScalarExp {
    x: PrimeFieldElem,
    next_exp_x: PrimeFieldElem,
}

impl Iterator for ScalarExp {
    type Item = PrimeFieldElem;

    fn next(&mut self) -> Option<PrimeFieldElem> {
        let exp_x = self.next_exp_x.clone();
        self.next_exp_x = &self.next_exp_x * &self.x;
        Some(exp_x)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (usize::max_value(), None)
    }
}

pub fn exp_iter_type1(x: &PrimeFieldElem) -> ScalarExp {
    let next_exp_x = x.f.elem(&1u8);
    ScalarExp { x: x.clone(), next_exp_x }
}

pub fn exp_iter_type2(x: &PrimeFieldElem) -> ScalarExp {
    let next_exp_x = x.clone();
    ScalarExp { x: x.clone(), next_exp_x }
}

pub fn scalar_exp_vartime(x: &PrimeFieldElem, mut n: u64) -> PrimeFieldElem {
    let mut result = x.f.elem(&1u8);
    let mut aux = x.clone(); // x, x^2, x^4, x^8, ...

    while n > 0 {
        let bit = n & 1;
        if bit == 1 {
            result = &result * &aux;
        }
        n = n >> 1;
        aux = &aux * &aux;
    }
    result
}

pub fn sum_of_powers_type1(x: &PrimeFieldElem, n: usize) -> PrimeFieldElem {
    if !n.is_power_of_two() {
        return sum_of_powers_slow_type1(x, n);
    }
    if n == 0 || n == 1 {
        return x.f.elem(&n);
    }
    let mut m = n;
    let mut result = x.f.elem(&1u8) + x;
    let mut factor = x.clone();

    while m > 2 {
        factor = &factor * &factor;
        result = &result + &factor * &result;
        m = m / 2;
    }
    result
}

fn sum_of_powers_slow_type1(x: &PrimeFieldElem, n: usize) -> PrimeFieldElem {
    let mut sum = x.f.elem(&0u8);
    for x in exp_iter_type1(x).take(n) {
        sum = sum + x;
    }
    sum
}

pub fn sum_of_powers_type2(x: &PrimeFieldElem, n: usize) -> PrimeFieldElem {
    if !n.is_power_of_two() {
        return sum_of_powers_slow_type2(x, n);
    }
    if n == 0 || n == 1 {
        return x.f.elem(&n);
    }
    let mut m = n;
    let mut result = x + x * x;
    let mut factor = x.clone();

    while m > 2 {
        factor = &factor * &factor;
        result = &result + &factor * &result;
        m = m / 2;
    }
    result
}

fn sum_of_powers_slow_type2(x: &PrimeFieldElem, n: usize) -> PrimeFieldElem {
    let mut sum = x.f.elem(&0u8);
    for x in exp_iter_type2(x).take(n) {
        sum = sum + x;
    }
    sum
}

// #[allow(dead_code)]
// pub fn inner_product(a: &[Scalar], b: &[Scalar]) -> Scalar {
//     let mut out = Scalar::zero();
//     for i in 0..a.len() {
//         out += a[i] * b[i];
//     }
//     out
// }

pub fn weighted_inner_product(
    curve: &Rc<Secp256k1>,
    a: &[PrimeFieldElem],
    b: &[PrimeFieldElem],
    c: &[PrimeFieldElem],
) -> PrimeFieldElem {
    let mut out = curve.f.elem(&0u8);
    for i in 0..a.len() {
        out = &out + (&a[i] * &b[i] * &c[i]);
    }
    out
}
