use crate::FieldExt;
use ecfft::{prepare_domain, prepare_matrices, GoodCurve, Matrix2x2};

#[derive(Clone, Debug)]
pub struct ECFFTConfig<F: FieldExt> {
    pub domain: Vec<Vec<F>>,
    pub matrices: Vec<Vec<Matrix2x2<F>>>,
    pub inverse_matrices: Vec<Vec<Matrix2x2<F>>>,
}

pub fn gen_config<F: FieldExt>(num_cols: usize) -> ECFFTConfig<F> {
    assert!(num_cols.is_power_of_two());
    let expansion_factor = 2;
    let codeword_len = num_cols * expansion_factor;

    let k = (codeword_len as f64).log2() as usize;

    let good_curve = GoodCurve::find_k(k);
    let domain = prepare_domain(good_curve);
    let (matrices, inverse_matrices) = prepare_matrices(&domain);

    ECFFTConfig {
        domain,
        matrices,
        inverse_matrices,
    }
}
