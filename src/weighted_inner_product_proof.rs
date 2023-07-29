#![allow(non_snake_case)]

extern crate alloc;

use alloc::vec::Vec;
use rand_core::OsRng;
use std::mem;

use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use curve25519_dalek::scalar::Scalar;
use merlin::Transcript;

use crate::util;
use crate::errors::ProofError;
use crate::publickey::PublicKey;
use crate::transcript::TranscriptProtocol;
use crate::mulvec::MulVec;

/**
 * Wieghted inner product proof
 * The size of the proof is
 *   2 * log_2{n} + 2 : CompressedRistretto,
 *   3 : Scalar
 */
#[derive(Clone, Debug)]
pub struct WeightedInnerProductProof {
    pub(crate) L_vec: Vec<CompressedRistretto>,
    pub(crate) R_vec: Vec<CompressedRistretto>,
    pub(crate) A: CompressedRistretto,
    pub(crate) B: CompressedRistretto,
    pub(crate) r_prime: Scalar,
    pub(crate) s_prime: Scalar,
    pub(crate) d_prime: Scalar,
}

impl WeightedInnerProductProof {
    /**
     * Prove weighted inner product
     */
    pub fn prove(
        transcript: &mut Transcript,
        pk: &PublicKey,
        a_vec: &Vec<Scalar>,
        b_vec: &Vec<Scalar>,
        power_of_y_vec: &Vec<Scalar>,
        gamma: Scalar,
        commitment: RistrettoPoint,
    ) -> Self {
        // random number generator
        let mut csprng = OsRng;

        // create slices G, H, a, b, c
        let mut G = &mut pk.G_vec.clone()[..];
        let mut H = &mut pk.H_vec.clone()[..];
        let mut a = &mut a_vec.clone()[..];
        let mut b = &mut b_vec.clone()[..];
        let mut power_of_y = &mut power_of_y_vec.clone()[..];

        // create copyed mutable scalars
        let mut alpha = gamma;

        // create copyed mutable commitment
        let mut P = commitment;

        // all of the input vectors must have the same length
        let mut n = G.len();
        assert_eq!(H.len(), n);
        assert_eq!(a.len(), n);
        assert_eq!(b.len(), n);
        assert_eq!(power_of_y.len(), n);

        // the length should be power of two
        assert!(n.is_power_of_two());

        // set transcript weight vector
        transcript.weighted_inner_product_domain_sep(power_of_y_vec);

        // allocate memory for L_vec and R_vec
        let logn = n.next_power_of_two().trailing_zeros() as usize;
        let mut L_vec: Vec<CompressedRistretto> = Vec::with_capacity(logn);
        let mut R_vec: Vec<CompressedRistretto> = Vec::with_capacity(logn);

        // n > 1 case
        while n != 1 {
            n = n / 2;

            // split a, b, c, G, H vector
            let (a1, a2) = a.split_at_mut(n);
            let (b1, b2) = b.split_at_mut(n);
            let (power_of_y1, power_of_y2) = power_of_y.split_at_mut(n);
            let (G1, G2) = G.split_at_mut(n);
            let (H1, H2) = H.split_at_mut(n);

            // compute c_L and c_R
            let c_L = util::weighted_inner_product(&a1, &b2, &power_of_y1);
            let c_R = util::weighted_inner_product(&a2, &b1, &power_of_y2);

            // random d_L and d_R by prover
            let d_L: Scalar = Scalar::random(&mut csprng);
            let d_R: Scalar = Scalar::random(&mut csprng);

            // compute L and R
            let y_nhat = power_of_y1[n - 1];
            let y_nhat_inv = y_nhat.invert();
            let G1_exp: Vec<Scalar> = a2.iter().map(|a2_i| y_nhat * a2_i).collect();
            let G2_exp: Vec<Scalar> = a1.iter().map(|a1_i| y_nhat_inv * a1_i).collect();

            let mut mv_g2 = MulVec::new();
            mv_g2.add_scalars(&G2_exp);
            mv_g2.add_scalars(&b2.to_vec());
            mv_g2.add_scalar(&c_L);
            mv_g2.add_scalar(&d_L);

            mv_g2.add_points(&G2.to_vec());
            mv_g2.add_points(&H1.to_vec());
            mv_g2.add_point(&pk.g);
            mv_g2.add_point(&pk.h);
            let L = mv_g2.calculate();

            let mut mv_g1 = MulVec::new();
            mv_g1.add_scalars(&G1_exp);
            mv_g1.add_scalars(&b1.to_vec());
            mv_g1.add_scalar(&c_R);
            mv_g1.add_scalar(&d_R);

            mv_g1.add_points(&G1.to_vec());
            mv_g1.add_points(&H2.to_vec());
            mv_g1.add_point(&pk.g);
            mv_g1.add_point(&pk.h);
            let R = mv_g1.calculate();

            L_vec.push(L.compress());
            R_vec.push(R.compress());

            // get challenge e
            transcript.append_point(b"L", &(L.compress()));
            transcript.append_point(b"R", &(R.compress()));
            let e = transcript.challenge_scalar(b"e");
            let e_inv = e.invert();
            let e_sqr = e * e;
            let e_sqr_inv = e_inv * e_inv;

            // update a, b, c, alpha, G, H, P
            let mut mv1 = MulVec::new();
            mv1.add_scalar(&e_sqr);
            mv1.add_scalar(&e_sqr_inv);
            mv1.add_point(&L);
            mv1.add_point(&R);
            mv1.calculate();
            P = P + mv1.calculate();

            let y_nhat_e_inv = y_nhat * e_inv;
            let y_nhat_inv_e = y_nhat_inv * e;

            for i in 0..n {
                a1[i] = a1[i] * e + a2[i] * y_nhat_e_inv;
                b1[i] = b1[i] * e_inv + b2[i] * e;

                let mut mv_g = MulVec::new();
                mv_g.add_scalar(&e_inv);
                mv_g.add_scalar(&y_nhat_inv_e);
                mv_g.add_point(&G1[i]);
                mv_g.add_point(&G2[i]);
                G1[i] = mv_g.calculate();

                let mut mv_h = MulVec::new();
                mv_h.add_scalar(&e);
                mv_h.add_scalar(&e_inv);
                mv_h.add_point(&H1[i]);
                mv_h.add_point(&H2[i]);
                H1[i] = mv_h.calculate();
            }
            a = a1;
            b = b1;
            power_of_y = power_of_y1;
            G = G1;
            H = H1;
            alpha += e_sqr * d_L + e_sqr_inv * d_R;
        }

        // random r, s, delta, eta
        let r: Scalar = Scalar::random(&mut csprng);
        let s: Scalar = Scalar::random(&mut csprng);
        let delta: Scalar = Scalar::random(&mut csprng);
        let eta: Scalar = Scalar::random(&mut csprng);

        // compute A and B
        let rcbsca = r * power_of_y[0] * b[0] + s * power_of_y[0] * a[0];
        let rcs = r * power_of_y[0] * s;

        let mut mv_A = MulVec::new();
        mv_A.add_scalar(&r);
        mv_A.add_scalar(&s);
        mv_A.add_scalar(&rcbsca);
        mv_A.add_scalar(&delta);
        mv_A.add_point(&G[0]);
        mv_A.add_point(&H[0]);
        mv_A.add_point(&pk.g);
        mv_A.add_point(&pk.h);
        let A = mv_A.calculate().compress();

        let mut mv_B = MulVec::new();
        mv_B.add_scalar(&rcs);
        mv_B.add_scalar(&eta);
        mv_B.add_point(&pk.g);
        mv_B.add_point(&pk.h);
        let B = mv_B.calculate().compress();

        // get challenge e
        transcript.append_point(b"A", &A);
        transcript.append_point(b"B", &B);
        let e = transcript.challenge_scalar(b"e");

        // compute r_prime, s_prime, delta_prime
        let r_prime = r + a[0] * e;
        let s_prime = s + b[0] * e;
        let d_prime = eta + delta * e + alpha * e * e;

        WeightedInnerProductProof {
            L_vec: L_vec,
            R_vec: R_vec,
            A: A,
            B: B,
            r_prime: r_prime,
            s_prime: s_prime,
            d_prime: d_prime,
        }
    }

    /**
     * To represent all verification process in one
     * multi-exponentiation, this function gets exponents
     * of commitment which can be computted publicly.
     *
     * Commitment = A' + Sum G_vec[i] * G_exp_of_commitment[i]
     *             + Sum H_vec[i] * H_exp_of_commitment[i]
     *             + g * g_exp_of_commitment + Sum V * V_exp_of_commitment
     */
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        pk: &PublicKey,
        power_of_y_vec: &Vec<Scalar>,
        G_exp_of_commitment: &[Scalar],
        H_exp_of_commitment: &[Scalar],
        g_exp_of_commitment: &Scalar,
        V_exp_of_commitment: &[Scalar],
        A_prime: RistrettoPoint,
        V: &[RistrettoPoint],
    ) -> Result<(), ProofError> {
        use curve25519_dalek::traits::IsIdentity;

        let logn = self.L_vec.len();
        let n = (1 << logn) as usize;
        let y = power_of_y_vec[0];

        let (challenges_sqr, challenges_inv_sqr, s_vec, e)
            = self.verification_scalars(n, power_of_y_vec, transcript)?;

        let s_prime_vec = s_vec.iter().rev();
        let e_sqr = e * e;
        let r_prime_e_y = self.r_prime * e * y;
        let s_prime_e = self.s_prime * e;

        // compute RHS / LHS
        let Ls_exp = challenges_sqr
            .iter()
            .map(|challenges_sqr_i| challenges_sqr_i * e_sqr)
            .collect();

        let Rs_exp = challenges_inv_sqr
            .iter()
            .map(|challenges_inv_sqr_i| challenges_inv_sqr_i * e_sqr)
            .collect();

        let G_exp = s_vec
            .iter()
            .zip(G_exp_of_commitment.iter())
            .zip(util::exp_iter_type2(y.invert()))
            .map(|((s_vec_i, g_exp_of_comm_i), power_of_y_vec_inv_i)| {
                -s_vec_i * power_of_y_vec_inv_i * r_prime_e_y + g_exp_of_comm_i * e_sqr
            })
            .collect();

        let H_exp = s_prime_vec.zip(H_exp_of_commitment.iter()).map(
            |(s_prime_vec_i, h_exp_of_comm_i)| {
                -s_prime_vec_i * s_prime_e + h_exp_of_comm_i * e_sqr
            },
        ).collect();

        let g_exp = -self.r_prime * y * self.s_prime + *g_exp_of_commitment * e_sqr;
        let h_exp = -self.d_prime;
        let V_exp = V_exp_of_commitment
            .iter()
            .map(|V_exp_of_commitment_i| V_exp_of_commitment_i * e_sqr)
            .collect();

        let mut mv = MulVec::new();
        mv.add_scalar(&Scalar::one());
        mv.add_scalar(&e);
        mv.add_scalar(&e_sqr);
        mv.add_scalar(&g_exp);
        mv.add_scalar(&h_exp);
        mv.add_scalars(&Ls_exp);
        mv.add_scalars(&Rs_exp);
        mv.add_scalars(&G_exp);
        mv.add_scalars(&H_exp);
        mv.add_scalars(&V_exp);

        let L_vec2 = self.L_vec.iter().map(|x| x.decompress().unwrap()).collect();
        let R_vec2 = self.R_vec.iter().map(|x| x.decompress().unwrap()).collect();

        mv.add_point(&self.B.decompress().unwrap());
        mv.add_point(&self.A.decompress().unwrap());
        mv.add_point(&A_prime);
        mv.add_point(&pk.g);
        mv.add_point(&pk.h);
        mv.add_points(&L_vec2);
        mv.add_points(&R_vec2);
        mv.add_points(&pk.G_vec);
        mv.add_points(&pk.H_vec);
        mv.add_points(&V.to_vec());

        let expected = mv.calculate();

        // check LSH == RHS
        if expected.is_identity() {
            Ok(())
        } else {
            Err(ProofError::VerificationError)
        }
    }

    pub fn batch_invert(xs: &Vec<Scalar>) -> (Scalar, Vec<Scalar>) {
        let mut prod = Scalar::one();
        let mut inv_xs = vec![];

        for x in xs {
            inv_xs.push(x.invert());
            prod = prod * x.invert();
        }

        (prod, inv_xs)
    }

    pub fn verification_scalars(
        &self,
        n: usize,
        power_of_y_vec: &Vec<Scalar>,
        transcript: &mut Transcript,
    ) -> Result<(Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Scalar), ProofError> {
        let logn = self.L_vec.len();
        if n != (1 << logn) {
            return Err(ProofError::VerificationError);
        }

        transcript.weighted_inner_product_domain_sep(power_of_y_vec);

        // 1. Recompute e_1, e_2, ..., e_k based on the proof transcript
        let mut challenges = Vec::with_capacity(logn);
        for (L, R) in self.L_vec.iter().zip(self.R_vec.iter()) {
            transcript.validate_and_append_point(b"L", L)?;
            transcript.validate_and_append_point(b"R", R)?;
            challenges.push(transcript.challenge_scalar(b"e"));
        }

        // 2. Compute 1/(e_1...e_k) and 1/e_k, ..., 1/e_1
        // let mut challenges_inv = challenges.clone();
        // let allinv = Scalar::batch_invert(&mut challenges_inv);
        let (allinv, mut challenges_inv) =
            WeightedInnerProductProof::batch_invert(&challenges);

        // 3. Compute e_i^2 and (1/e_i)^2
        for i in 0..logn {
            challenges[i] = challenges[i] * challenges[i];
            challenges_inv[i] = challenges_inv[i] * challenges_inv[i];
        }
        let challenges_sqr = challenges;
        let challenges_inv_sqr = challenges_inv;

        // 4. Recompute e
        transcript.validate_and_append_point(b"A", &self.A)?;
        transcript.validate_and_append_point(b"B", &self.B)?;
        let e = transcript.challenge_scalar(b"e");

        // 5. Compute s_vec and s_prime_vec
        let mut s_vec = Vec::with_capacity(n);
        s_vec.push(allinv);

        for i in 1..n {
            let log_i = (32 - 1 - (i as u32).leading_zeros()) as usize;
            let k = 1 << log_i;
            let u_log_i_sq = challenges_sqr[(logn - 1) - log_i];
            s_vec.push(s_vec[i - k] * u_log_i_sq);
        }

        Ok((challenges_sqr, challenges_inv_sqr, s_vec, e))
    }
    //
    pub fn size(&self) -> usize {
        let n = self.L_vec.len();
        let mut res: usize = 0;
        for i in 0..n {
            res += mem::size_of_val(self.L_vec[i].as_bytes());
            res += mem::size_of_val(self.R_vec[i].as_bytes());
        }
        res += mem::size_of_val(self.A.as_bytes());
        res += mem::size_of_val(self.B.as_bytes());
        res += mem::size_of_val(self.r_prime.as_bytes());
        res += mem::size_of_val(self.s_prime.as_bytes());
        res += mem::size_of_val(self.d_prime.as_bytes());
        res
    }
}
