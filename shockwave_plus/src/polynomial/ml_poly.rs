use tensor_pcs::EqPoly;

use crate::FieldExt;

#[derive(Clone, Debug)]
pub struct MlPoly<F> {
    pub evals: Vec<F>,
    pub num_vars: usize,
}

impl<F: FieldExt> MlPoly<F> {
    pub fn new(evals: Vec<F>) -> Self {
        assert!(evals.len().is_power_of_two());
        let num_vars = (evals.len() as f64).log2() as usize;
        Self { evals, num_vars }
    }

    fn dot_prod(x: &[F], y: &[F]) -> F {
        assert_eq!(x.len(), y.len());
        let mut result = F::ZERO;
        for i in 0..x.len() {
            result += x[i] * y[i];
        }
        result
    }

    // Evaluate the multilinear extension of the polynomial `a`, at point `t`.
    // `a` is in evaluation form.
    pub fn eval(&self, t: &[F]) -> F {
        let n = self.evals.len();
        debug_assert_eq!((n as f64).log2() as usize, t.len());
        // Evaluate the multilinear extension of the polynomial `a`,
        // over the boolean hypercube

        let eq_evals = EqPoly::new(t.to_vec()).evals();

        Self::dot_prod(&self.evals, &eq_evals)
    }
}
