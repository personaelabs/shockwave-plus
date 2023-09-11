use crate::polynomial::eq_poly::EqPoly;
use crate::r1cs::r1cs::Matrix;
use crate::sumcheck::SumCheckProof;
use crate::tensor_pcs::hasher::Hasher;
use crate::tensor_pcs::TensorMultilinearPCS;
use crate::transcript::TranscriptLike;
use crate::FieldGC;

use super::sumcheck::prove_sum;

pub struct SumCheckPhase2<F: FieldGC> {
    A_mat: Matrix<F>,
    B_mat: Matrix<F>,
    C_mat: Matrix<F>,
    Z_evals: Vec<F>,
    rx: Vec<F>,
    r: [F; 3],
}

impl<F: FieldGC> SumCheckPhase2<F> {
    pub fn new(
        A_mat: Matrix<F>,
        B_mat: Matrix<F>,
        C_mat: Matrix<F>,
        Z_evals: Vec<F>,
        rx: Vec<F>,
        r: [F; 3],
    ) -> Self {
        Self {
            A_mat,
            B_mat,
            C_mat,
            Z_evals,
            rx,
            r,
        }
    }

    pub fn prove<H: Hasher<F>>(
        &self,
        pcs: &TensorMultilinearPCS<F, H>,
        transcript: &mut impl TranscriptLike<F>,
        blind: bool,
    ) -> SumCheckProof<F, H> {
        let r_A = self.r[0];
        let r_B = self.r[1];
        let r_C = self.r[2];

        let n = self.Z_evals.len();
        let num_vars = (self.Z_evals.len() as f64).log2() as usize;

        let evals_rx = EqPoly::new(self.rx.clone()).evals();
        let mut A_evals = vec![F::ZERO; n];
        let mut B_evals = vec![F::ZERO; n];
        let mut C_evals = vec![F::ZERO; n];

        for entry in &self.A_mat.entries {
            A_evals[entry.col] += evals_rx[entry.row] * entry.val;
        }
        for entry in &self.B_mat.entries {
            B_evals[entry.col] += evals_rx[entry.row] * entry.val;
        }
        for entry in &self.C_mat.entries {
            C_evals[entry.col] += evals_rx[entry.row] * entry.val;
        }

        let mut eval_tables = vec![
            A_evals.clone(),
            B_evals.clone(),
            C_evals.clone(),
            self.Z_evals.clone(),
        ];

        let poly_degree = 2;
        let comb_func = |x: &[F]| (x[0] * r_A + x[1] * r_B + x[2] * r_C) * x[3];

        prove_sum(
            num_vars,
            poly_degree,
            &mut eval_tables,
            comb_func,
            blind,
            pcs,
            transcript,
            "sc_phase_2".to_string(),
        )
    }
}
