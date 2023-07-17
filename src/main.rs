use curve25519_dalek::scalar::Scalar;
use merlin::Transcript;
use bulletproofsplus::PublicKey;
use bulletproofsplus::range::{RangeProof, RangeProver};

fn main() {
    let n = 64;
    let m = 2;
    let pk = PublicKey::new(n * m);
    let mut prover = RangeProver::new();

    let v1 = 2u64;
    let v2 = 5u64;
    let gamma1 = Scalar::from(3u8);
    let gamma2 = Scalar::from(7u8);
    prover.commit(&pk, v1, gamma1);
    prover.commit(&pk, v2, gamma2);

    let mut prover_transcript = Transcript::new(b"RangeProof Test");
    let proof: RangeProof = RangeProof::prove(
        &mut prover_transcript,
        &pk,
        n,
        &prover,
    );
    let commitment_vec = prover.commitment_vec;
    let mut transcript = Transcript::new(b"RangeProof Test");
    let result = proof.verify(
        &mut transcript,
        &pk,
        n,
        &commitment_vec,
    );
    assert_eq!(result, Ok(()));
    println!("worked");
}