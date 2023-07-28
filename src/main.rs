use bulletproofsplus::PublicKey;
#[allow(unused_imports)]
use bulletproofsplus::range::{RangeProof, RangeProver};
use bulletproofsplus::secp256k1::building_block::secp256k1::secp256k1::Secp256k1;
use std::rc::Rc;

fn main() {
    println!("started");

    let n = 64;
    let m = 1;

    let curve = Rc::new(Secp256k1::new());
    println!("created secp256k1 curve");

    // setup generators
    let pk = PublicKey::new(&curve, n * m);
    println!("created public key");

    // create m commitmentments and add to prover
    let mut prover = RangeProver::new(curve.clone());
    let v1 = 2u64;
    let gamma1 = curve.f_n.elem(&3u8);
    prover.commit(&pk, v1, gamma1);

    // let v2 = 5u64;
    // let gamma2 = curve.f_n.elem(&7u8);
    // prover.commit(&pk, v2, gamma2);

    // build proof
    println!("started proving...");
    let proof: RangeProof = RangeProof::prove(
        curve.clone(),
        &pk,
        n,
        &prover,
    );
    println!("finished proving.");

    let commitment_vec = prover.commitment_vec;

    // verify proof
    println!("started verifying...");
    let result = proof.verify(
        curve.clone(),
        &pk,
        n,
        &commitment_vec,
    );
    println!("finished verifying.");

    assert_eq!(result, Ok(()));
}