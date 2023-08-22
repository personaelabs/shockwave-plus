#![allow(non_snake_case)]
use criterion::{criterion_group, criterion_main, Criterion};
use shockwave_plus::good_curves::secp256k1::secp256k1_good_curve;
use shockwave_plus::{det_num_cols, ShockwavePlus, Transcript, R1CS};

fn shockwave_plus_bench(c: &mut Criterion) {
    type F = halo2curves::secp256k1::Fp;

    for exp in [12, 15, 18] {
        let num_cons = 2usize.pow(exp as u32);
        let num_input = 3;
        let num_vars = num_cons - num_input;

        let (r1cs, witness) = R1CS::<F>::produce_synthetic_r1cs(num_vars, num_input);

        let mut group = c.benchmark_group(format!("ShockwavePlus num_cons: {}", r1cs.num_cons()));
        let l = 319;
        let num_cols = det_num_cols(r1cs.z_len(), l);

        let (good_curve, coset_offset) =
            secp256k1_good_curve((num_cols as f64).log2() as usize + 1);

        group.bench_function("config", |b| {
            b.iter(|| {
                ShockwavePlus::new(r1cs.clone(), l, good_curve, coset_offset);
            })
        });

        let shockwave_plus = ShockwavePlus::new(r1cs.clone(), l, good_curve, coset_offset);

        group.bench_function("prove", |b| {
            b.iter(|| {
                let mut transcript = Transcript::new(b"bench");
                shockwave_plus.prove(&witness, &r1cs.public_input, &mut transcript);
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
