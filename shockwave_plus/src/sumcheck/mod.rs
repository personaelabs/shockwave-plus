mod sc_phase_1;
mod sc_phase_2;
use crate::{FieldGC, TensorMLOpening};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
pub use sc_phase_1::SumCheckPhase1;
pub use sc_phase_2::SumCheckPhase2;
pub mod unipoly;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SumCheckProof<F: FieldGC> {
    pub round_poly_coeffs: Vec<Vec<F>>,
    pub blinder_poly_sum: Option<F>,
    pub blinder_poly_eval_proof: Option<TensorMLOpening<F>>,
}

impl<F: FieldGC> SumCheckProof<F> {
    pub fn is_blinded(&self) -> bool {
        self.blinder_poly_sum.is_some()
    }
}
