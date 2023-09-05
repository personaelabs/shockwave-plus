use crate::polynomial::eq_poly::EqPoly;

use crate::FieldGC;

#[derive(Clone, Debug)]
pub struct MlPoly<F> {
    pub evals: Vec<F>,
    pub num_vars: usize,
}

impl<F: FieldGC> MlPoly<F> {
    #[allow(dead_code)]
    pub fn new(evals: Vec<F>) -> Self {
        assert!(evals.len().is_power_of_two());
        let num_vars = (evals.len() as f64).log2() as usize;
        Self { evals, num_vars }
    }

    #[allow(dead_code)]
    fn dot_prod(x: &[F], y: &[F]) -> F {
        assert_eq!(x.len(), y.len());
        let mut result = F::from(0u64);
        for i in 0..x.len() {
            result += x[i] * y[i];
        }
        result
    }

    // Evaluate the multilinear extension of the polynomial `a`, at point `t`.
    // `a` is in evaluation form.
    // `t` should be in big-endian form.
    #[allow(dead_code)]
    pub fn eval(&self, t: &[F]) -> F {
        let n = self.evals.len();
        debug_assert_eq!((n as f64).log2() as usize, t.len());

        let eq_evals = EqPoly::new(t.to_vec()).evals();

        Self::dot_prod(&self.evals, &eq_evals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = ark_secp256k1::Fq;

    #[test]
    fn test_ml_poly_eval() {
        let ZERO = F::from(0u32);
        let ONE = F::from(1u32);

        let num_vars = 4;
        let num_evals = 2usize.pow(num_vars as u32);
        let evals = (0..num_evals)
            .map(|x| F::from(x as u64))
            .collect::<Vec<F>>();

        let ml_poly = MlPoly::new(evals.clone());
        let eval_last = ml_poly.eval(&[ONE, ONE, ONE, ONE]);
        assert_eq!(
            eval_last,
            evals[evals.len() - 1],
            "The last evaluation is not correct"
        );

        let eval_first = ml_poly.eval(&[ZERO, ZERO, ZERO, ZERO]);
        assert_eq!(eval_first, evals[0], "The first evaluation is not correct");

        let eval_second = ml_poly.eval(&[ZERO, ZERO, ZERO, ONE]);
        assert_eq!(
            eval_second, evals[1],
            "The second evaluation is not correct"
        );
    }
}
