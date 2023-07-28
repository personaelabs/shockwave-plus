use crate::FieldExt;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct UniPoly<F: FieldExt> {
    pub coeffs: Vec<F>,
}

impl<F: FieldExt> UniPoly<F> {
    fn eval_cubic(&self, x: F) -> F {
        // ax^3 + bx^2 + cx + d
        let x_sq = x.square();
        let x_cub = x_sq * x;

        let a = self.coeffs[0];
        let b = self.coeffs[1];
        let c = self.coeffs[2];
        let d = self.coeffs[3];

        a * x_cub + b * x_sq + c * x + d
    }

    fn eval_quadratic(&self, x: F) -> F {
        // ax^3 + bx^2 + cx + d
        let x_sq = x.square();

        let a = self.coeffs[0];
        let b = self.coeffs[1];
        let c = self.coeffs[2];

        a * x_sq + b * x + c
    }

    pub fn eval(&self, x: F) -> F {
        if self.coeffs.len() == 3 {
            self.eval_quadratic(x)
        } else {
            self.eval_cubic(x)
        }
    }

    pub fn interpolate(evals: &[F]) -> Self {
        debug_assert!(
            evals.len() == 4 || evals.len() == 3,
            "Only cubic and quadratic polynomials are supported"
        );

        let two_inv = F::TWO_INV;

        if evals.len() == 4 {
            // ax^3 + bx^2 + cx + d
            let six_inv = F::from(6u64).invert().unwrap();

            let d = evals[0];
            let a = six_inv
                * (evals[3] - evals[2] - evals[2] - evals[2] + evals[1] + evals[1] + evals[1]
                    - evals[0]);
            let b = two_inv
                * (evals[0] + evals[0] - evals[1] - evals[1] - evals[1] - evals[1] - evals[1]
                    + evals[2]
                    + evals[2]
                    + evals[2]
                    + evals[2]
                    - evals[3]);

            let c = evals[1] - d - a - b;

            Self {
                coeffs: vec![a, b, c, d],
            }
        } else {
            let c = evals[0];
            let a = (evals[2] - evals[1] - evals[1] + evals[0]) * two_inv;
            let b = evals[1] - a - c;

            Self {
                coeffs: vec![a, b, c],
            }
        }
    }
}
