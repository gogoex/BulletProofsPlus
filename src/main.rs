use rand_core::OsRng;
use curve25519_dalek::scalar::Scalar;
use merlin::Transcript;
use bulletproofsplus::PublicKey;
use bulletproofsplus::range::{RangeProof, RangeProver};

fn main() {
    let n = 8;
    let m = 1;
    let pk = PublicKey::new(n * m);
    let mut prover = RangeProver::new();

    for _i in 0..m {
        prover.commit(&pk, 31u64, Scalar::random(&mut OsRng));
    }

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

    println!("succeeded");
}