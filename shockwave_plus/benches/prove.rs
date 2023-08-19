#![allow(non_snake_case)]
use criterion::{criterion_group, criterion_main, Criterion};
use shockwave_plus::ShockwavePlus;
use shockwave_plus::R1CS;
use tensor_pcs::rs_config;
use tensor_pcs::rs_config::good_curves::secp256k1::secp256k1_good_curve;
use tensor_pcs::TensorRSMultilinearPCSConfig;
use tensor_pcs::{det_num_cols, Transcript};

fn shockwave_plus_bench(c: &mut Criterion) {
    type F = halo2curves::secp256k1::Fp;

    for exp in [12, 15, 18] {
        let num_cons = 2usize.pow(exp as u32);
        let num_input = 3;
        let num_vars = num_cons - num_input;

        let (r1cs, witness, pub_input) = R1CS::<F>::produce_synthetic_r1cs(num_vars, num_input);

        let mut group = c.benchmark_group(format!("ShockwavePlus num_cons: {}", r1cs.num_cons()));
        let l = 319;
        let num_cols = det_num_cols(r1cs.z_len(), l);

        let (good_curve, coset_offset) =
            secp256k1_good_curve((num_cols as f64).log2() as usize + 1);

        group.bench_function("config", |b| {
            b.iter(|| {
                rs_config::ecfft::gen_config_form_curve(good_curve, coset_offset);
            })
        });

        let ecfft_config = rs_config::ecfft::gen_config_form_curve(good_curve, coset_offset);
        let pcs_config = TensorRSMultilinearPCSConfig {
            expansion_factor: 2,
            l,
            ecfft_config: Some(ecfft_config),
            fft_domain: None,
            domain_powers: None,
        };

        let shockwave_plus = ShockwavePlus::new(r1cs.clone(), pcs_config);

        group.bench_function("prove", |b| {
            b.iter(|| {
                let mut transcript = Transcript::new(b"bench");
                shockwave_plus.prove(&witness, &pub_input, &mut transcript);
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
