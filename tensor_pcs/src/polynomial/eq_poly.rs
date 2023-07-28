use crate::FieldExt;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct EqPoly<F: FieldExt> {
    t: Vec<F>,
}

impl<F: FieldExt> EqPoly<F> {
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
        let ell = self.t.len(); // 4

        let mut evals: Vec<F> = vec![F::ONE; 2usize.pow(ell as u32)];
        let mut size = 1;
        for j in 0..ell {
            // in each iteration, we double the size of chis
            size *= 2; // 2 4 8 16
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2]; // i = 0, 2, 4, 7
                evals[i] = scalar * self.t[j]; // (1 * t0)(1 * t1)
                evals[i - 1] = scalar - evals[i]; // 1 - (1 * t0)(1 * t1)
            }
        }
        evals
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::sparse_ml_poly::SparseMLPoly;
    use halo2curves::ff::Field;

    type F = halo2curves::secp256k1::Fp;

    pub fn dot_prod<F: FieldExt>(x: &[F], y: &[F]) -> F {
        assert_eq!(x.len(), y.len());
        let mut result = F::ZERO;
        for i in 0..x.len() {
            result += x[i] * y[i];
        }
        result
    }

    #[test]
    fn test_eq_poly() {
        let m = 4;
        let t = (0..m).map(|i| F::from((i + 33) as u64)).collect::<Vec<F>>();
        let eq_poly = EqPoly::new(t.clone());
        eq_poly.evals();
    }
}
