#![allow(non_snake_case)]
use criterion::{criterion_group, criterion_main, Criterion};
use shockwave_plus::ShockwavePlus;
use shockwave_plus::R1CS;
use tensor_pcs::Transcript;

fn shockwave_plus_bench(c: &mut Criterion) {
    type F = halo2curves::secp256k1::Fp;

    for exp in [12, 15, 18] {
        let num_vars = 2usize.pow(exp);
        let num_input = 3;

        let (r1cs, witness) = R1CS::<F>::produce_synthetic_r1cs(num_vars, num_input);

        let mut group = c.benchmark_group(format!("ShockwavePlus num_cons: {}", r1cs.num_cons));
        let l = 319;
        let num_rows = (((2f64 / l as f64).sqrt() * (num_vars as f64).sqrt()) as usize)
            .next_power_of_two()
            / 2;
        let ShockwavePlus = ShockwavePlus::new(r1cs.clone(), l, num_rows);
        group.bench_function("prove", |b| {
            b.iter(|| {
                let mut transcript = Transcript::new(b"bench");
                ShockwavePlus.prove(&witness, &r1cs.public_input, &mut transcript);
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
