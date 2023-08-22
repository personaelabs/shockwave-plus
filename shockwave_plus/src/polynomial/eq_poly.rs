use crate::FieldExt;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct EqPoly<F: FieldExt> {
    t: Vec<F>,
}

impl<F: FieldExt> EqPoly<F> {
    // `t` should be in big-endian.
    pub fn new(t: Vec<F>) -> Self {
        Self { t }
    }

    pub fn eval(&self, x: &[F]) -> F {
        let mut result = F::ONE;
        let one = F::ONE;

        for i in 0..x.len() {
            result *= self.t[i] * x[i] + (one - self.t[i]) * (one - x[i]);
        }
        result
    }

    // Copied from microsoft/Spartan
    pub fn evals(&self) -> Vec<F> {
        let ell = self.t.len();

        let mut evals: Vec<F> = vec![F::ONE; 2usize.pow(ell as u32)];
        let mut size = 1;
        for j in 0..ell {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2];
                evals[i] = scalar * self.t[j];
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = halo2curves::secp256k1::Fp;
    use halo2curves::ff::Field;

    #[test]
    fn test_eq_poly() {
        let m = 4;
        let t = (0..m).map(|i| F::from((i + 33) as u64)).collect::<Vec<F>>();
        let eq_poly = EqPoly::new(t.clone());
        let evals = eq_poly.evals();

        let eval_first = eq_poly.eval(&[F::ZERO, F::ZERO, F::ZERO, F::ZERO]);
        assert_eq!(eval_first, evals[0], "The first evaluation is not correct");

        let eval_second = eq_poly.eval(&[F::ZERO, F::ZERO, F::ZERO, F::ONE]);
        assert_eq!(
            eval_second, evals[1],
            "The second evaluation is not correct"
        );

        let eval_last = eq_poly.eval(&[F::ONE, F::ONE, F::ONE, F::ONE]);
        assert_eq!(
            eval_last,
            evals[evals.len() - 1],
            "The last evaluation is not correct"
        );
    }
}
