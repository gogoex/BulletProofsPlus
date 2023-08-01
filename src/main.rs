use bulletproofsplus::{PublicKey, bls12_381::building_block::scalar::prime_field_elem::PrimeFieldElem};
#[allow(unused_imports)]
use bulletproofsplus::range::{RangeProof, RangeProver};
use bulletproofsplus::bls12_381::building_block::arith::Arith;

fn main() {
    println!("started");
    Arith::init();

    let n = 64;
    let m = 2;

    // setup generators
    let pk = PublicKey::new(n * m);
    println!("created public key");

    // create m commitmentments and add to prover
    let mut prover = RangeProver::new();
    let v1 = 2u64;
    let gamma1 = PrimeFieldElem::new(3);
    prover.commit(&pk, v1, gamma1);

    let v2 = 5u64;
    let gamma2 = PrimeFieldElem::new(7);
    prover.commit(&pk, v2, gamma2);

    // build proof
    println!("started proving...");
    let proof: RangeProof = RangeProof::prove(
        &pk,
        n,
        &prover,
    );
    println!("finished proving.");

// println!("Ls={:?}", &proof.proof.L_vec);
// println!("Rs={:?}", &proof.proof.R_vec);
// println!("A={:?}", &proof.A);
// println!("A_wip={:?}", &proof.proof.A);
// println!("B={:?}", &proof.proof.B);
// println!("r'={:?}", &proof.proof.r_prime);
// println!("s'={:?}", &proof.proof.s_prime);
// println!("delta'={:?}", &proof.proof.d_prime);

    let commitment_vec = prover.commitment_vec;

    // verify proof
    println!("started verifying...");
    let result = proof.verify(
        &pk,
        n,
        &commitment_vec,
    );
    println!("finished verifying.");

    assert_eq!(result, Ok(()));
}