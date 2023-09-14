#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::EqPoly;
use crate::FieldGC;

#[derive(Clone, Debug)]
pub struct SparseMLPoly<F> {
    pub evals: Vec<(usize, F)>,
    pub num_vars: usize,
}

impl<F: FieldGC> SparseMLPoly<F> {
    pub fn new(evals: Vec<(usize, F)>, num_vars: usize) -> Self {
        Self { evals, num_vars }
    }

    pub fn num_entries(&self) -> usize {
        2usize.pow(self.num_vars as u32)
    }

    // `t` should be in big-endian form.
    pub fn eval(&self, t: &[F]) -> F {
        debug_assert_eq!(self.num_vars, t.len());
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

    pub fn eval_naive(&self, t: &[F]) -> F {
        debug_assert_eq!(self.num_vars, t.len());
        // Evaluate the multilinear extension of the polynomial `a`,
        // over the boolean hypercube

        let eq_poly = EqPoly::new(t.to_vec());

        #[cfg(feature = "parallel")]
        let result = self
            .evals
            .par_iter()
            .map(|eval| eq_poly.eval_as_bits(eval.0) * eval.1)
            .sum();

        #[cfg(not(feature = "parallel"))]
        let result = self
            .evals
            .iter()
            .map(|eval| eq_poly.eval_as_bits(eval.0) * eval.1)
            .sum();

        result
    }
}
