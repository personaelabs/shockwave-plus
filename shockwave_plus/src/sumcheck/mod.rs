mod sc_phase_1;
mod sc_phase_2;
use crate::{tensor_pcs::hasher::Hasher, FieldGC, TensorMLOpening};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
pub use sc_phase_1::SumCheckPhase1;
pub use sc_phase_2::SumCheckPhase2;
pub mod sumcheck;
pub mod unipoly;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct SumCheckProof<F: FieldGC, H: Hasher<F>> {
    pub label: String,
    pub round_poly_coeffs: Vec<Vec<F>>,
    pub blinder_poly_sum: Option<F>,
    pub blinder_poly_eval_proof: Option<TensorMLOpening<F, H>>,
}

impl<F: FieldGC, H: Hasher<F>> SumCheckProof<F, H> {
    pub fn is_blinded(&self) -> bool {
        self.blinder_poly_sum.is_some()
    }
}
