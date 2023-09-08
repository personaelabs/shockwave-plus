#![allow(non_snake_case)]
use criterion::{criterion_group, criterion_main, Criterion};
use shockwave_plus::{
    IOPattern, PoseidonCurve, PoseidonHasher, PoseidonTranscript, TensorRSMultilinearPCSConfig,
};
use shockwave_plus::{ShockwavePlus, R1CS};

fn shockwave_plus_bench(c: &mut Criterion) {
    type F = ark_secp256k1::Fq;

    for exp in [12, 15, 18] {
        let num_cons = 2usize.pow(exp as u32);
        let num_input = 3;
        let num_vars = num_cons - num_input;

        let (r1cs, witness, pub_input) = R1CS::<F>::produce_synthetic_r1cs(num_vars, num_input);

        let mut group = c.benchmark_group(format!("ShockwavePlus num_cons: {}", r1cs.num_cons()));

        let l = 319;
        let expansion_factor = 2;
        let blind = true;
        let pcs_config = TensorRSMultilinearPCSConfig::new(r1cs.z_len(), expansion_factor, l);
        let poseidon_hasher = PoseidonHasher::new(PoseidonCurve::SECP256K1);
        let shockwave_plus = ShockwavePlus::new(r1cs.clone(), pcs_config, poseidon_hasher);

        group.bench_function("prove", |b| {
            b.iter(|| {
                let mut transcript = PoseidonTranscript::new(
                    b"bench",
                    PoseidonCurve::SECP256K1,
                    IOPattern::new(vec![]),
                );
                shockwave_plus.prove(&witness, &pub_input, &mut transcript, blind);
            })
        });

        let proof = shockwave_plus
            .prove(
                &witness,
                &pub_input,
                &mut PoseidonTranscript::new(
                    b"bench",
                    PoseidonCurve::SECP256K1,
                    IOPattern::new(vec![]),
                ),
                blind,
            )
            .0;

        group.bench_function("verify", |b| {
            b.iter(|| {
                let mut transcript = PoseidonTranscript::new(
                    b"bench",
                    PoseidonCurve::SECP256K1,
                    IOPattern::new(vec![]),
                );
                shockwave_plus.verify(&proof, &mut transcript);
            })
        });
    }
}

fn set_duration() -> Criterion {
    Criterion::default().sample_size(10)
}

criterion_group! {
    name = benches;
    config = set_duration();
    targets = shockwave_plus_bench
}
criterion_main!(benches);
