#![allow(non_snake_case)]

extern crate alloc;

use alloc::vec::Vec;

use crate::{
    errors::ProofError,
    publickey::PublicKey,
    bls12_381::building_block::{
        point::point::Point,
        scalar::prime_field_elem::PrimeFieldElem,
        mulvec::MulVec,
        zero::Zero2,
    },
    util,
    weighted_inner_product_proof::WeightedInnerProductProof,
};

pub mod prover;

pub use self::prover::RangeProver;

#[derive(Clone, Debug)]
pub struct RangeProof {
    pub A: Point,
    pub proof: WeightedInnerProductProof,
}

impl RangeProof {
    pub fn prove(
        pk: &PublicKey,
        n: usize,
        prover: &RangeProver,
    ) -> Self {
        let m = prover.v_vec.len();
        if m == 1 {
            // single proof case
            Self::prove_single(
                pk,
                n,
                prover.v_vec[0],
                &prover.gamma_vec[0],
                &prover.commitment_vec[0])
        } else {
            Self::prove_multiple(
                pk,
                n,
                m,
                &prover.v_vec,
                &prover.gamma_vec,
                &prover.commitment_vec,
            )
        }
    }

    pub fn verify(
        &self,
        pk: &PublicKey,
        n: usize,
        commitment_vec: &[Point],
    ) -> Result<(), ProofError> {
        let m = commitment_vec.len();
        if m == 1 {
            self.verify_single(
                pk,
                n,
                &commitment_vec[0],
            )
        } else {
            self.verify_multiple(
                pk,
                n,
                m,
                commitment_vec,
            )
        }
    }

    fn prove_single(
        pk: &PublicKey,
        n: usize,
        v: u64,
        gamma: &PrimeFieldElem,
        commitment: &Point,
    ) -> RangeProof {
        println!("in prove single");

        // check parameter
        assert_eq!(pk.G_vec.len(), n);
        assert_eq!(pk.H_vec.len(), n);

        // random alpha
        let alpha = &PrimeFieldElem::new(7);

        // compute A
        let mut v_bits: Vec<u8> = Vec::with_capacity(n);
        let mut A = &pk.h * alpha;
        let mut i = 0;
        for (G_i, H_i) in pk.G_vec.iter().zip(pk.H_vec.iter()) {
            v_bits.push(((v >> i) & 1) as u8);
            let mut point = -H_i.clone();
            if v_bits[i] != 0 { point = G_i.clone() };
            A = A + point;
            i += 1;
        }

        // get challenges
        let y = &PrimeFieldElem::new(7);
        let z = &PrimeFieldElem::new(7);

        // compute A_hat
        let one = &PrimeFieldElem::new(1);
        let two = &PrimeFieldElem::new(2);
        let power_of_two: Vec<PrimeFieldElem> = util::exp_iter_type1(&PrimeFieldElem::new(2)).take(n).collect();
        let power_of_y: Vec<PrimeFieldElem> = util::exp_iter_type2(y).take(n).collect();
        let power_of_y_rev = power_of_y.clone().into_iter().rev();

        let mut G_vec_sum = Point::zero();
        for p in &pk.G_vec {
            G_vec_sum = G_vec_sum + p;
        }

        let G_vec_sum_exp = -z;
        let H_exp: Vec<PrimeFieldElem> = power_of_two
            .iter()
            .zip(power_of_y_rev)
            .map(|(power_of_two_i, power_of_y_rev_i)| power_of_two_i * power_of_y_rev_i + z)
            .collect();

        let V_exp = util::scalar_exp_vartime(&y, (n + 1) as u64);

        let mut g_exp = PrimeFieldElem::new(0);
        for x in &power_of_y {
            g_exp = g_exp + x;
        }
        g_exp = g_exp * (z - (z * z));
        g_exp = g_exp - ((util::scalar_exp_vartime(&two, n as u64) - one) * &V_exp * z);

        let mut mv = MulVec::new();
        mv.add_scalar(&PrimeFieldElem::new(1));
        mv.add_scalar(&G_vec_sum_exp);
        mv.add_scalars(&H_exp);
        mv.add_scalar(&g_exp);
        mv.add_scalar(&V_exp);

        mv.add_point(&A);
        mv.add_point(&G_vec_sum);
        mv.add_points(&pk.H_vec);
        mv.add_point(&pk.g);
        mv.add_point(&commitment);

        let A_hat = mv.calculate();

        // compute a_vec, b_vec, alpha_hat
        let nz = -z;
        let one_minus_z = one - z;

        let a_vec: Vec<PrimeFieldElem> = v_bits
            .iter()
            .map(|v_bits_i| if *v_bits_i == 0 { nz.clone() } else { one_minus_z.clone() })
            .collect();

        let b_vec: Vec<PrimeFieldElem> = H_exp
            .iter()
            .zip(v_bits.iter())
            .map(|(H_exp_i, v_bits_i)| {
                if *v_bits_i == 0 { (H_exp_i - one).clone() } else { H_exp_i.clone() }
            })
            .collect();

        let alpha_hat = alpha + gamma * &V_exp;

        // generate weighted inner product proof
        let proof = WeightedInnerProductProof::prove(
            &pk,
            &a_vec,
            &b_vec,
            &power_of_y,
            alpha_hat,
            A_hat,
        );
        RangeProof {
            A: A,
            proof: proof,
        }
    }

    fn verify_single(
        &self,
        pk: &PublicKey,
        n: usize,
        commitment: &Point,
    ) -> Result<(), ProofError> {
        println!("in verify single");

        // get challenges
        let y = &PrimeFieldElem::new(7);
        let z = &PrimeFieldElem::new(7);

        // decompress A
        let As = self.A.clone();
        let Vs = commitment.clone();

        // compute exponent of A_hat
        let one = &PrimeFieldElem::new(1);
        let two = &PrimeFieldElem::new(2);

        let power_of_two: Vec<PrimeFieldElem> = util::exp_iter_type1(&PrimeFieldElem::new(2)).take(n).collect();
        let power_of_y: Vec<PrimeFieldElem> = util::exp_iter_type2(y).take(n).collect();
        let power_of_y_rev = power_of_y.iter().rev();

        let G_exp: Vec<PrimeFieldElem> = vec![-z; n];
        let H_exp: Vec<PrimeFieldElem> = power_of_two
            .iter()
            .zip(power_of_y_rev)
            .map(|(power_of_two_i, power_of_y_rev_i)| power_of_two_i * power_of_y_rev_i + z)
            .collect();
        let V_exp = util::scalar_exp_vartime(&y, (n + 1) as u64);

        let mut g_exp = PrimeFieldElem::new(0);
        for x in &power_of_y {
            g_exp = g_exp + x;
        }
        g_exp = g_exp * (z - (z * z));
        g_exp = g_exp - ((util::scalar_exp_vartime(&two, n as u64) - one) * &V_exp * z);

        self.proof.verify(
            &pk,
            &power_of_y,
            &G_exp,
            &H_exp,
            &g_exp,
            &[V_exp],
            As,
            &[Vs]
        )
    }

    fn prove_multiple(
        pk: &PublicKey,
        n: usize,
        m: usize,
        v: &[u64],
        gamma_vec: &[PrimeFieldElem],
        commitment_vec: &[Point],
    ) -> RangeProof {
        println!("in prove multiple");

        let mn = n * m;
        // check parameter
        assert_eq!(pk.G_vec.len(), mn);
        assert_eq!(pk.H_vec.len(), mn);

        // random alpha
        let alpha = &PrimeFieldElem::new(33);

        // compute A
        let mut v_bits: Vec<u8> = Vec::with_capacity(mn);
        let mut A = &pk.h * alpha;
        let mut i = 0;

        for (G_i, H_i) in pk.G_vec.iter().zip(pk.H_vec.iter()) {
            let index1 = i % n;
            let index2 = i / n;

            v_bits.push(((v[index2] >> index1) & 1) as u8);

            let mut point = -H_i;

            if v_bits[i] != 0 {
                point = G_i.clone();
            }

            A = A + point;
            i += 1;
        }
        let y = &PrimeFieldElem::new(12);
        let z = &PrimeFieldElem::new(23);
// println!("aL={:?}", &v_bits);

        // compute d
        let power_of_two: Vec<PrimeFieldElem> = util::exp_iter_type1(&PrimeFieldElem::new(2)).take(n).collect();
        let power_of_y: Vec<PrimeFieldElem> = util::exp_iter_type2(&y).take(mn).collect();
        let power_of_y_rev = power_of_y.iter().rev();
        let z_sqr = &(z * z);
        let power_of_z: Vec<PrimeFieldElem> = util::exp_iter_type2(&z_sqr).take(m).collect();
        let d: Vec<PrimeFieldElem> = power_of_z
            .iter()
            .flat_map(|exp_z| power_of_two.iter().map(move |exp_2| exp_2 * exp_z))
            .collect();
// println!("power_of_z={:?}", &power_of_z);
// println!("d={:?}", &d);

        // compute A_hat
        let G_vec_sum_exp = &-z;

        let H_exp: Vec<PrimeFieldElem> = d
            .iter()
            .zip(power_of_y_rev)
            .map(|(d_i, power_of_y_rev_i)| d_i * power_of_y_rev_i + z)
            .collect();

        let power_of_y_mn_plus_1 = &util::scalar_exp_vartime(&y, (mn + 1) as u64);

        let V_exp: Vec<PrimeFieldElem> = power_of_z
            .iter()
            .map(|power_of_z_i| power_of_z_i * power_of_y_mn_plus_1)
            .collect();

        let mut g_exp = PrimeFieldElem::new(0);
        for x in &power_of_y {
            g_exp = &g_exp + x;
        }
        g_exp = g_exp * (z - z_sqr);

        let mut d_sum = PrimeFieldElem::new(0);
        for x in d {
            d_sum = &d_sum + x;
        }
        g_exp = &g_exp - (&d_sum * power_of_y_mn_plus_1 * z);

        // let G_vec_sum: AffinePoint = pk.G_vec.iter().sum();
        let mut G_vec_sum = Point::zero();
        for x in &pk.G_vec {
            G_vec_sum = &G_vec_sum + x;
        }
//println!("A={:?}", &(&A * PrimeFieldElem::new(1)));

        let mut mv = MulVec::new();
        mv.add_scalar(&PrimeFieldElem::new(1));
        mv.add_scalar(&G_vec_sum_exp);
        mv.add_scalars(&H_exp);
        mv.add_scalar(&g_exp);
        mv.add_scalars(&V_exp);

        mv.add_point(&A);
        mv.add_point(&G_vec_sum);
        mv.add_points(&pk.H_vec);
        mv.add_point(&pk.g);
        mv.add_points(&commitment_vec.to_vec());

        let A_hat = mv.calculate();
// println!("A_hat={:?}", &A_hat);

        // compute a_vec, b_vec, alpha_hat
        let one = &PrimeFieldElem::new(1);
        let nz = &-z;
        let one_minus_z = &(one - z);

        let a_vec: Vec<PrimeFieldElem> = v_bits
            .iter()
            //.map(|v_bits_i| Scalar::conditional_select(&nz, &one_minus_z, *v_bits_i))
            .map(|v_bits_i| if *v_bits_i != 0 { one_minus_z.clone() } else { nz.clone() })
            .collect();

        let b_vec: Vec<PrimeFieldElem> = H_exp
            .iter()
            .zip(v_bits.iter())
            .map(|(H_exp_i, v_bits_i)| {
                //Scalar::conditional_select(&(H_exp_i - one), H_exp_i, *v_bits_i)
                if *v_bits_i != 0 { H_exp_i.clone() } else { (H_exp_i - one).clone() }
            })
            .collect();

        let power_of_z_gamma_vec: Vec<PrimeFieldElem> = power_of_z
            .iter()
            .zip(gamma_vec.iter())
            .map(|(power_of_z_i, gamma_i)| power_of_z_i * gamma_i )
            .collect();
        let mut power_of_z_gamma_sum = PrimeFieldElem::new(0);
        for x in power_of_z_gamma_vec {
            power_of_z_gamma_sum = power_of_z_gamma_sum + x;
        }

        let alpha_hat = alpha + power_of_z_gamma_sum * power_of_y_mn_plus_1;
// println!("alpha_hat={:?}", &alpha_hat);

// println!("a_hat_R={:?}", &b_vec);

        // generate weighted inner product proof
        let proof = WeightedInnerProductProof::prove(
            &pk,
            &a_vec,
            &b_vec,
            &power_of_y,
            alpha_hat,
            A_hat,
        );

// println!("Ls={:?}", &proof.L_vec);
// println!("Rs={:?}", &proof.R_vec);
//println!("A_wip={:?}", &proof.A);
// println!("B={:?}", &proof.B);
// println!("r_prime={:?}", &proof.r_prime);
// println!("s_prime={:?}", &proof.s_prime);
// println!("d_prime={:?}", &proof.d_prime);

        RangeProof {
            A,
            proof,
        }
    }

    fn verify_multiple(
        &self,
        pk: &PublicKey,
        n: usize,
        m: usize,
        commitment_vec: &[Point],
    ) -> Result<(), ProofError> {
        println!("in verify multiple");

        let mn = n * m;

        // 1. Recompute y and z
        let y = &PrimeFieldElem::new(12);
        let z = &PrimeFieldElem::new(23);
        let minus_z = &-z;
        let z_sqr = &(z * z);

        // 2. Compute power of two, power of y, power of z
        let power_of_two: Vec<PrimeFieldElem> = util::exp_iter_type1(&PrimeFieldElem::new(2)).take(n).collect();
        let mut power_of_y: Vec<PrimeFieldElem> = util::exp_iter_type2(&y).take(mn + 1).collect();
        let power_of_y_mn_plus_1 = match power_of_y.pop() {
            Some(point) => point,
            None => panic!("fail to pop"),
        };
        let mut power_of_y_rev: Vec<PrimeFieldElem> = power_of_y.clone();
        power_of_y_rev.reverse();

        let power_of_z: Vec<PrimeFieldElem> = util::exp_iter_type2(&z_sqr).take(m).collect();

        // 3. Compute concat_z_and_2
        let concat_z_and_2: Vec<PrimeFieldElem> = power_of_z
            .iter()
            .flat_map(|exp_z| power_of_two.iter().map(move |exp_2| exp_2 * exp_z))
            .collect();

        // 4. Compute scalars for verification
        let (challenges_sqr, challenges_inv_sqr, s_vec, e)
            = self.proof.verification_scalars(mn)?;

        let s_prime_vec = s_vec.iter().rev();
        let e_inv = &e.inv();
        let e_sqr = &(&e * &e);
        let e_sqr_inv = &e_sqr.inv();
        let r_prime_e_inv_y = &(&self.proof.r_prime * e_inv * y);
        let s_prime_e_inv = &(&self.proof.s_prime * e_inv);

        // 5. Compute exponents of G_vec, H_vec, g, and h
        let r_prime = &self.proof.r_prime;
        let s_prime = &self.proof.s_prime;
        let d_prime = &self.proof.d_prime;

        let G_exp: Vec<PrimeFieldElem> = s_vec.iter()
            .zip(util::exp_iter_type2(&y.inv()))
            .map(|(s_vec_i, power_of_y_inv_i)| minus_z - s_vec_i * power_of_y_inv_i * r_prime_e_inv_y)
            .collect();

        let H_exp: Vec<PrimeFieldElem> = s_prime_vec
            .zip(concat_z_and_2.iter())
            .zip(power_of_y_rev)
            .map(|((s_prime_vec_i, d_i), power_of_y_rev_i)| - s_prime_e_inv * s_prime_vec_i + (d_i * power_of_y_rev_i + z))
            .collect();

        let sum_y = util::sum_of_powers_type2(&y, mn);
        let sum_2 = util::sum_of_powers_type1(&PrimeFieldElem::new(2), n);
        let sum_z = util::sum_of_powers_type2(&z_sqr, m);

        let g_exp = -r_prime * s_prime * y * e_sqr_inv + (sum_y * (z - z_sqr) - &power_of_y_mn_plus_1 * z * sum_2 * sum_z);
        let h_exp = -d_prime * e_sqr_inv;

        // 6. Compute exponents of V_vec
        let V_exp: Vec<PrimeFieldElem> = power_of_z.iter()
            .map(|power_of_z_i| power_of_z_i * &power_of_y_mn_plus_1)
            .collect();

        // 7. Compute RHS / LHS
        let mut mv = MulVec::new();
        mv.add_scalar(&PrimeFieldElem::new(1));
        mv.add_scalar(&e_inv);
        mv.add_scalar(&e_sqr_inv);
        mv.add_scalar(&g_exp);
        mv.add_scalar(&h_exp);
        mv.add_scalars(&challenges_sqr);
        mv.add_scalars(&challenges_inv_sqr);
        mv.add_scalars(&G_exp);
        mv.add_scalars(&H_exp);
        mv.add_scalars(&V_exp);

        mv.add_point(&self.A);
        mv.add_point(&self.proof.A);
        mv.add_point(&self.proof.B);
        mv.add_point(&pk.g);
        mv.add_point(&pk.h);
        mv.add_points(&self.proof.L_vec);
        mv.add_points(&self.proof.R_vec);
        mv.add_points(&pk.G_vec);
        mv.add_points(&pk.H_vec);
        mv.add_points(&commitment_vec.to_vec());

        let expected = mv.calculate();

        if expected.is_zero() {
            Ok(())
        } else {
            Err(ProofError::VerificationError)
        }
    }

    // pub fn size(&self) -> usize {
    //     let mut res: usize = 0;
    //     res += mem::size_of_val(self.A.as_bytes());
    //     res += self.proof.size();
    //     res
    // }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     fn range_proof(
//         n: usize,
//         m: usize,
//     ) {
//         let pk = PublicKey::new(n * m);
//         let curve = Rc::new(Secp256k1::new());
//         let mut prover = RangeProver::new(curve);

//         for _i in 0..m {
//             prover.commit(&pk, 31u64, curve.f_n.rand_elem(true));
//         }

//         let proof: RangeProof = RangeProof::prove(
//             curve.clone(),
//             &pk,
//             n,
//             &prover,
//         );

//         let mut verifier = RangeVerifier::new();
//         verifier.allocate(&prover.commitment_vec);

//         let result = proof.verify(
//             curve.clone(),
//             &pk,
//             n,
//             &verifier.commitment_vec,
//         );
//         assert_eq!(result, Ok(()));
//     }

//     #[test]
//     fn test_range_proof_all() {
//         range_proof(32 as usize, 1 as usize);
//         range_proof(32 as usize, 2 as usize);
//         range_proof(32 as usize, 4 as usize);
//         range_proof(32 as usize, 8 as usize);
//         range_proof(64 as usize, 1 as usize);
//         range_proof(64 as usize, 2 as usize);
//         range_proof(64 as usize, 4 as usize);
//         range_proof(64 as usize, 8 as usize);
//     }
// }