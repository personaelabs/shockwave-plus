use crate::FieldExt;
use ecfft::{prepare_domain, prepare_matrices, GoodCurve, Matrix2x2};

#[derive(Clone, Debug)]
pub struct ECFFTConfig<F: FieldExt> {
    pub domain: Vec<Vec<F>>,
    pub matrices: Vec<Vec<Matrix2x2<F>>>,
    pub inverse_matrices: Vec<Vec<Matrix2x2<F>>>,
}

pub fn gen_config_form_curve<F: FieldExt>(
    good_curve: GoodCurve<F>,
    coset_offset: (F, F),
) -> ECFFTConfig<F> {
    let domain = prepare_domain(good_curve, coset_offset.0, coset_offset.1);
    let (matrices, inverse_matrices) = prepare_matrices(&domain);

    ECFFTConfig {
        domain,
        matrices,
        inverse_matrices,
    }
}
