use crate::FieldExt;

pub struct BlinderPoly<F: FieldExt> {
    inner_poly_coeffs: Vec<Vec<F>>,
}

impl<F: FieldExt> BlinderPoly<F> {
    pub fn sample_random(num_vars: usize, degree: usize) -> Self {
        let mut rng = rand::thread_rng();
        let inner_poly_coeffs = (0..num_vars)
            .map(|_| (0..(degree + 1)).map(|_| F::random(&mut rng)).collect())
            .collect();

        Self { inner_poly_coeffs }
    }

    pub fn eval(&self, x: &[F]) -> F {
        let mut res = F::ZERO;

        for (coeffs, x_i) in self.inner_poly_coeffs.iter().zip(x.iter()) {
            let mut tmp = F::ZERO;
            let mut x_i_pow = F::ONE;

            for coeff in coeffs.iter() {
                tmp += *coeff * x_i_pow;
                x_i_pow *= x_i;
            }

            res += tmp;
        }

        res
    }
}
