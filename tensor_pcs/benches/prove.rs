use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tensor_pcs::{
    rs_config, FieldExt, SparseMLPoly, TensorMultilinearPCS, TensorRSMultilinearPCSConfig,
    Transcript,
};

fn poly<F: FieldExt>(num_vars: usize) -> SparseMLPoly<F> {
    let num_entries: usize = 2usize.pow(num_vars as u32);

    let evals = (0..num_entries)
        .map(|i| (i, F::from(i as u64)))
        .collect::<Vec<(usize, F)>>();

    let ml_poly = SparseMLPoly::new(evals, num_vars);
    ml_poly
}

fn config_base<F: FieldExt>(ml_poly: &SparseMLPoly<F>) -> TensorRSMultilinearPCSConfig<F> {
    let num_vars = ml_poly.num_vars;
    let num_evals = 2usize.pow(num_vars as u32);
    let num_rows = 2usize.pow((num_vars / 2) as u32);

    let expansion_factor = 2;

    TensorRSMultilinearPCSConfig::<F> {
        expansion_factor,
        domain_powers: None,
        fft_domain: None,
        ecfft_config: None,
        l: 10,
        num_entries: num_evals,
        num_rows,
    }
}

fn pcs_fft_bench(c: &mut Criterion) {
    type F = halo2curves::pasta::Fp;

    let num_vars = 13;
    let ml_poly = poly(num_vars);
    let open_at = (0..ml_poly.num_vars)
        .map(|i| F::from(i as u64))
        .collect::<Vec<F>>();

    let mut config = config_base(&ml_poly);
    config.fft_domain = Some(rs_config::smooth::gen_config::<F>(config.num_cols()));

    let mut group = c.benchmark_group("pcs fft");
    group.bench_function("prove", |b| {
        b.iter(|| {
            let pcs = TensorMultilinearPCS::<F>::new(config.clone());

            let mut transcript = Transcript::new(b"bench");
            let comm = pcs.commit(&black_box(ml_poly.clone()));
            pcs.open(&comm, &ml_poly, &open_at, &mut transcript);
        })
    });
}

fn pcs_ecfft_bench(c: &mut Criterion) {
    type F = halo2curves::secp256k1::Fp;

    let num_vars = 13;
    let ml_poly = poly(num_vars);
    let open_at = (0..ml_poly.num_vars)
        .map(|i| F::from(i as u64))
        .collect::<Vec<F>>();

    let mut config = config_base(&ml_poly);
    config.ecfft_config = Some(rs_config::ecfft::gen_config::<F>(config.num_cols()));

    let mut group = c.benchmark_group("pcs ecfft");
    group.bench_function("prove", |b| {
        b.iter(|| {
            let pcs = TensorMultilinearPCS::<F>::new(config.clone());

            let mut transcript = Transcript::new(b"bench");
            let comm = pcs.commit(&black_box(ml_poly.clone()));
            pcs.open(&comm, &ml_poly, &open_at, &mut transcript);
        })
    });
}

fn set_duration() -> Criterion {
    Criterion::default().sample_size(10)
}

criterion_group! {
    name = benches;
    config = set_duration();
    targets = pcs_fft_bench, pcs_ecfft_bench
}
criterion_main!(benches);
