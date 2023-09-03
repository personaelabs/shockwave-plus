use crate::EqPoly;
use ark_ff::PrimeField;

#[derive(Clone, Debug)]
pub struct SparseMLPoly<F> {
    pub evals: Vec<(usize, F)>,
    pub num_vars: usize,
}

impl<F: PrimeField> SparseMLPoly<F> {
    pub fn new(evals: Vec<(usize, F)>, num_vars: usize) -> Self {
        Self { evals, num_vars }
    }

    pub fn num_entries(&self) -> usize {
        2usize.pow(self.num_vars as u32)
    }

    pub fn from_dense(dense_evals: Vec<F>) -> Self {
        assert!(dense_evals.len().is_power_of_two());
        let sparse_evals = dense_evals
            .iter()
            .enumerate()
            .filter(|(_, eval)| **eval != F::ZERO)
            .map(|(i, eval)| (i, *eval))
            .collect::<Vec<(usize, F)>>();
        let num_vars = (dense_evals.len() as f64).log2() as usize;

        Self {
            evals: sparse_evals,
            num_vars,
        }
    }

    // `t` should be in big-endian form.
    pub fn eval(&self, t: &[F]) -> F {
        assert_eq!(self.num_vars, t.len());
        // Evaluate the multilinear extension of the polynomial `a`,
        // over the boolean hypercube

        let eq_poly = EqPoly::new(t.to_vec());
        let eq_evals = eq_poly.evals();

        let mut result = F::ZERO;

        for eval in &self.evals {
            result += eq_evals[eval.0] * eval.1;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = ark_secp256k1::Fq;

    #[test]
    fn test_sparse_ml_poly_eval() {
        let ZERO = F::from(0u32);
        let ONE = F::from(1u32);

        let num_vars = 4;
        let num_evals = 2usize.pow(num_vars as u32);
        let evals = (0..num_evals)
            .map(|x| F::from((x as u64) as u64))
            .collect::<Vec<F>>();

        let ml_poly = SparseMLPoly::from_dense(evals.clone());
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
