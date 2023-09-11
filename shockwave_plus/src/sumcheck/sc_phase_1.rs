use crate::polynomial::eq_poly::EqPoly;
use crate::sumcheck::SumCheckProof;
use crate::tensor_pcs::hasher::Hasher;
use crate::tensor_pcs::TensorMultilinearPCS;
use crate::transcript::TranscriptLike;

use crate::FieldGC;

use super::sumcheck::prove_sum;

pub struct SumCheckPhase1<F: FieldGC> {
    Az_evals: Vec<F>,
    Bz_evals: Vec<F>,
    Cz_evals: Vec<F>,
}

impl<F: FieldGC> SumCheckPhase1<F> {
    pub fn new(Az_evals: Vec<F>, Bz_evals: Vec<F>, Cz_evals: Vec<F>) -> Self {
        Self {
            Az_evals,
            Bz_evals,
            Cz_evals,
        }
    }

    pub fn prove<H: Hasher<F>>(
        &self,
        pcs: &TensorMultilinearPCS<F, H>,
        transcript: &mut impl TranscriptLike<F>,
        blind: bool,
    ) -> (SumCheckProof<F, H>, (F, F, F)) {
        let poly_num_vars = (self.Az_evals.len() as f64).log2() as usize;
        let poly_degree = 3;

        let tau = transcript.challenge_vec(poly_num_vars, "tau".to_string());

        let mut eval_tables = vec![
            self.Az_evals.clone(),
            self.Bz_evals.clone(),
            self.Cz_evals.clone(),
            EqPoly::new(tau).evals(),
        ];
        let comb_func = |x: &[F]| (x[0] * x[1] - x[2]) * x[3];

        let sumcheck_proof = prove_sum(
            poly_num_vars,
            poly_degree,
            &mut eval_tables,
            comb_func,
            blind,
            pcs,
            transcript,
            "sc_phase_1".to_string(),
        );

        let v_A = eval_tables[0][0];
        let v_B = eval_tables[1][0];
        let v_C = eval_tables[2][0];

        (sumcheck_proof, (v_A, v_B, v_C))
    }
}
