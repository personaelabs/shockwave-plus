use crate::{EqPoly, FieldExt};

#[derive(Clone, Debug)]
pub struct SparseMLPoly<F> {
    pub evals: Vec<(usize, F)>,
    pub num_vars: usize,
}

impl<F: FieldExt> SparseMLPoly<F> {
    pub fn new(evals: Vec<(usize, F)>, num_vars: usize) -> Self {
        Self { evals, num_vars }
    }

    pub fn from_dense(dense_evals: Vec<F>) -> Self {
        let sparse_evals = dense_evals
            .iter()
            .filter(|eval| **eval != F::ZERO)
            .enumerate()
            .map(|(i, eval)| (i, *eval))
            .collect::<Vec<(usize, F)>>();
        let num_vars = (dense_evals.len() as f64).log2() as usize;

        Self {
            evals: sparse_evals,
            num_vars,
        }
    }

    pub fn eval(&self, t: &[F]) -> F {
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
