#![allow(non_snake_case)]
mod polynomial;
mod poseidon;
mod r1cs;
mod sumcheck;
mod tensor_pcs;
mod transcript;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer};
use polynomial::{eq_poly::EqPoly, ml_poly::MlPoly, sparse_ml_poly::SparseMLPoly};
use sumcheck::{SumCheckPhase1, SumCheckPhase2, SumCheckProof};

// Exports
pub use ark_ff;
pub use ark_secp256k1;
pub use ark_serialize;
pub use ecfft;
pub use poseidon::sponge::*;
pub use poseidon::{constants as poseidon_constants, Poseidon, PoseidonConstants, PoseidonCurve};
pub use r1cs::{Matrix, SparseMatrixEntry, R1CS};
pub use rs_config::good_curves::FieldGC;
pub use tensor_pcs::hasher::{Hasher, PoseidonHasher};
pub use tensor_pcs::rs_config::good_curves;
pub use tensor_pcs::{
    det_num_cols, det_num_rows, ecfft::GoodCurve, rs_config, TensorMLOpening, TensorMultilinearPCS,
    TensorRSMultilinearPCSConfig,
};
pub use transcript::{AppendToTranscript, PoseidonTranscript, TranscriptLike};

use crate::sumcheck::sumcheck::verify_sum;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: FieldGC, H: Hasher<F>> {
    pub pub_input: Vec<F>,
    pub z_comm: H::T,
    pub sc_proof_1: SumCheckProof<F, H>,
    pub sc_proof_2: SumCheckProof<F, H>,
    pub z_eval_proof: TensorMLOpening<F, H>,
    pub v_A: F,
    pub v_B: F,
    pub v_C: F,
}

impl<F: FieldGC, H: Hasher<F>> Proof<F, H> {
    pub fn is_blinded(&self) -> bool {
        self.sc_proof_1.is_blinded()
    }
}

#[derive(Clone)]
pub struct ShockwavePlus<F: FieldGC, H: Hasher<F>> {
    r1cs: R1CS<F>,
    pcs: TensorMultilinearPCS<F, H>,
}

impl<F: FieldGC, H: Hasher<F>> ShockwavePlus<F, H> {
    pub fn new(r1cs: R1CS<F>, config: TensorRSMultilinearPCSConfig<F>, hasher: H) -> Self {
        let ecfft_config = &config.ecfft_config;
        let curve_k = (ecfft_config.domain[0].len() as f64).log2() as usize;

        let min_num_entries = r1cs.num_vars.next_power_of_two();
        let min_num_cols = config.num_cols(min_num_entries);

        let max_num_entries = r1cs.z_len().next_power_of_two();
        let max_num_cols = config.num_cols(max_num_entries);
        // Make sure that there are enough columns to run the l queries
        assert!(min_num_cols > config.l);

        // Make sure that the FFTree is large enough
        assert!(curve_k > (max_num_cols as f64).log2() as usize);

        let pcs = TensorMultilinearPCS::new(config, hasher);

        Self { r1cs, pcs }
    }

    pub fn prove(
        &self,
        r1cs_witness: &[F],
        r1cs_input: &[F],
        transcript: &mut impl TranscriptLike<F>,
        blind: bool,
    ) -> (Proof<F, H>, Vec<F>) {
        // Multilinear extension requires the number of evaluations
        // to be a power of two to uniquely determine the polynomial
        let mut padded_r1cs_witness = r1cs_witness.to_vec();
        padded_r1cs_witness.resize(padded_r1cs_witness.len().next_power_of_two(), F::ZERO);
        let witness_poly = MlPoly::new(padded_r1cs_witness.clone());

        let Z = R1CS::construct_z(r1cs_witness, r1cs_input);

        // Commit the witness polynomial
        let comm_witness_timer = start_timer!(|| "Commit witness");
        let committed_witness = self.pcs.commit(&padded_r1cs_witness, blind);
        let witness_comm = committed_witness.committed_tree.root;
        end_timer!(comm_witness_timer);

        // Add the witness commitment to the transcript
        witness_comm.append_to_transcript(transcript);

        // ############################
        // Phase 1
        // ###################

        let m = (self.r1cs.z_len() as f64).log2() as usize;

        let mut Az_poly = self.r1cs.A.mul_vector(&Z);
        let mut Bz_poly = self.r1cs.B.mul_vector(&Z);
        let mut Cz_poly = self.r1cs.C.mul_vector(&Z);

        Az_poly.resize(Z.len(), F::ZERO);
        Bz_poly.resize(Z.len(), F::ZERO);
        Cz_poly.resize(Z.len(), F::ZERO);

        // Prove that the
        // Q(t) = \sum_{x \in {0, 1}^m} (Az_poly(x) * Bz_poly(x) - Cz_poly(x)) eq(t, x)
        // is a zero-polynomial using the sum-check protocol.
        // We evaluate Q(t) at $\tau$ and check that it is zero.

        let sc_phase_1_timer = start_timer!(|| "Sumcheck phase 1");

        let sc_phase_1 = SumCheckPhase1::new(Az_poly.clone(), Bz_poly.clone(), Cz_poly.clone());
        let (sc_proof_1, (v_A, v_B, v_C)) = sc_phase_1.prove(&self.pcs, transcript, blind);

        end_timer!(sc_phase_1_timer);

        transcript.append_fe(v_A);
        transcript.append_fe(v_B);
        transcript.append_fe(v_C);

        // Phase 2
        let r = transcript.challenge_vec(3, "r".to_string());

        // T_2 should equal teh evaluations of the random linear combined polynomials

        let rx = (0..m)
            .map(|i| transcript.get(&format!("sc_phase_1-challenge-{}", i)))
            .collect::<Vec<F>>();
        let sc_phase_2_timer = start_timer!(|| "Sumcheck phase 2");
        let sc_phase_2 = SumCheckPhase2::new(
            self.r1cs.A.clone(),
            self.r1cs.B.clone(),
            self.r1cs.C.clone(),
            Z.clone(),
            rx.clone(),
            r.as_slice().try_into().unwrap(),
        );

        let sc_proof_2 = sc_phase_2.prove(&self.pcs, transcript, blind);

        let ry = (0..m)
            .map(|i| transcript.get(&format!("sc_phase_2-challenge-{}", i)))
            .collect::<Vec<F>>();

        end_timer!(sc_phase_2_timer);

        let z_open_timer = start_timer!(|| "Open witness poly");
        // Prove the evaluation of the polynomial Z(y) at ry
        let z_eval_proof = self.pcs.open(
            &committed_witness,
            &padded_r1cs_witness,
            &ry[1..],
            witness_poly.eval(&ry[1..]),
            transcript,
            blind,
        );
        end_timer!(z_open_timer);

        // Prove the evaluation of the polynomials A(y), B(y), C(y) at ry

        let rx_ry = vec![ry, rx].concat();
        (
            Proof {
                pub_input: r1cs_input.to_vec(),
                z_comm: witness_comm,
                sc_proof_1,
                sc_proof_2,
                z_eval_proof,
                v_A,
                v_B,
                v_C,
            },
            rx_ry,
        )
    }

    pub fn verify(&self, proof: &Proof<F, H>, transcript: &mut impl TranscriptLike<F>) {
        proof.z_comm.append_to_transcript(transcript);

        let mle_timer = start_timer!(|| "ML Extension");
        let A_mle = self.r1cs.A.to_ml_extension();
        let B_mle = self.r1cs.B.to_ml_extension();
        let C_mle = self.r1cs.C.to_ml_extension();
        end_timer!(mle_timer);

        let m = (self.r1cs.z_len() as f64).log2() as usize;

        // ############################
        // Verify phase 1 sumcheck
        // ############################

        let tau = transcript.challenge_vec(m, "tau".to_string());

        let sc_phase1_sum_target = F::ZERO;

        // The final eval should equal
        let v_A = proof.v_A;
        let v_B = proof.v_B;
        let v_C = proof.v_C;

        let T_1_eq = EqPoly::new(tau);

        let sc_phase1_poly = |challenge: &[F]| (v_A * v_B - v_C) * T_1_eq.eval(challenge);

        verify_sum(
            &proof.sc_proof_1,
            &self.pcs,
            sc_phase1_sum_target,
            sc_phase1_poly,
            transcript,
        );

        // ############################
        // Verify phase 2 sumcheck
        // ############################

        transcript.append_fe(v_A);
        transcript.append_fe(v_B);
        transcript.append_fe(v_C);

        let r = transcript.challenge_vec(3, "r".to_string());
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let sc_phase2_sum_target = r_A * v_A + r_B * v_B + r_C * v_C;

        let rx = (0..m)
            .map(|i| transcript.get(&format!("sc_phase_1-challenge-{}", i)))
            .collect::<Vec<F>>();

        let sc_phase2_poly = |ry: &[F]| {
            let rx_ry = [&rx, ry].concat();
            let witness_eval = proof.z_eval_proof.y;

            let eval_timer = start_timer!(|| "Eval R1CS");
            let A_eval = A_mle.eval_naive(&rx_ry);
            let B_eval = B_mle.eval_naive(&rx_ry);
            let C_eval = C_mle.eval_naive(&rx_ry);
            end_timer!(eval_timer);

            let input = (0..self.r1cs.num_input)
                .map(|i| (i + 1, proof.pub_input[i]))
                .collect::<Vec<(usize, F)>>();

            let input_poly =
                SparseMLPoly::new(vec![vec![(0, F::ONE)], input].concat(), ry.len() - 1);
            let input_poly_eval = input_poly.eval(&ry[1..]);

            let z_eval = (F::ONE - ry[0]) * input_poly_eval + ry[0] * witness_eval;

            (r_A * A_eval + r_B * B_eval + r_C * C_eval) * z_eval
        };

        verify_sum(
            &proof.sc_proof_2,
            &self.pcs,
            sc_phase2_sum_target,
            sc_phase2_poly,
            transcript,
        );

        let pcs_verify_timer = start_timer!(|| "Verify PCS");
        self.pcs.verify(&proof.z_eval_proof, transcript);
        end_timer!(pcs_verify_timer);
    }
}

#[cfg(test)]
mod tests {

    use crate::{tensor_pcs::hasher::PoseidonHasher, transcript::PoseidonTranscript};

    use super::*;

    #[test]
    fn test_shockwave_plus() {
        type F = ark_secp256k1::Fq;

        let num_vars = 20;
        let num_input = 3;
        let l = 2;
        let expansion_factor = 2;

        let (r1cs, witness, pub_input) = R1CS::<F>::produce_synthetic_r1cs(num_vars, num_input);

        let config = TensorRSMultilinearPCSConfig::new(r1cs.z_len(), expansion_factor, l);

        let poseidon_hasher = PoseidonHasher::new(PoseidonCurve::SECP256K1);
        // Prove and verify with and without zero-knowledge
        let shockwave_plus = ShockwavePlus::new(r1cs.clone(), config, poseidon_hasher);

        for blind in [true, false] {
            let mut prover_transcript =
                PoseidonTranscript::new(b"test", PoseidonCurve::SECP256K1, IOPattern::new(vec![]));

            let (proof, _) =
                shockwave_plus.prove(&witness, &pub_input, &mut prover_transcript, blind);

            let mut verifier_transcript =
                PoseidonTranscript::new(b"test", PoseidonCurve::SECP256K1, IOPattern::new(vec![]));
            shockwave_plus.verify(&proof, &mut verifier_transcript);
        }
    }
}
