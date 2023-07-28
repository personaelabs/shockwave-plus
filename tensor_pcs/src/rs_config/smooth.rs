use crate::FieldExt;

pub fn gen_config<F: FieldExt>(num_cols: usize) -> Vec<F> {
    assert!(num_cols.is_power_of_two());
    let expansion_factor = 2;
    let codeword_len = num_cols * expansion_factor;

    let domain_generator = F::ROOT_OF_UNITY.pow(&[
        2u32.pow(32 - ((codeword_len as f64).log2() as u32)) as u64,
        0,
        0,
        0,
    ]);

    // Compute the FFT domain
    let mut fft_domain = Vec::with_capacity(codeword_len);
    fft_domain.push(F::ONE);
    for i in 0..(codeword_len - 1) {
        fft_domain.push(fft_domain[i] * domain_generator);
    }

    fft_domain
}
