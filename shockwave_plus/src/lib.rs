#![allow(non_snake_case)]
mod polynomial;
#[cfg(feature = "poseidon")]
mod poseidon;
mod r1cs;
mod sumcheck;
mod tensor_pcs;
mod transcript;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer};
use polynomial::{eq_poly::EqPoly, ml_poly::MlPoly, sparse_ml_poly::SparseMLPoly};
use sumcheck::{SCPhase1Proof, SCPhase2Proof, SumCheckPhase1, SumCheckPhase2};

// Exports
pub use ark_ff;
pub use ark_secp256k1;
pub use ark_serialize;
#[cfg(feature = "poseidon")]
pub use poseidon::{constants as poseidon_constants, Poseidon, PoseidonConstants, PoseidonCurve};
pub use r1cs::{Matrix, SparseMatrixEntry, R1CS};
pub use rs_config::good_curves::FieldGC;
pub use tensor_pcs::rs_config::good_curves;
pub use tensor_pcs::{
    det_num_cols, det_num_rows, ecfft::GoodCurve, rs_config, TensorMLOpening, TensorMultilinearPCS,
    TensorRSMultilinearPCSConfig,
};
pub use transcript::{AppendToTranscript, Transcript};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSpartanProof<F: FieldGC> {
    pub pub_input: Vec<F>,
    pub z_comm: [u8; 32],
    pub sc_proof_1: SCPhase1Proof<F>,
    pub sc_proof_2: SCPhase2Proof<F>,
    pub z_eval_proof: TensorMLOpening<F>,
    pub v_A: F,
    pub v_B: F,
    pub v_C: F,
}

pub struct FullSpartanProof<F: FieldGC> {
    pub partial_proof: PartialSpartanProof<F>,
    pub A_eval_proof: TensorMLOpening<F>,
    pub B_eval_proof: TensorMLOpening<F>,
    pub C_eval_proof: TensorMLOpening<F>,
}

#[derive(Clone)]
pub struct ShockwavePlus<F: FieldGC> {
    r1cs: R1CS<F>,
    pcs: TensorMultilinearPCS<F>,
}

impl<F: FieldGC> ShockwavePlus<F> {
    pub fn new(r1cs: R1CS<F>, config: TensorRSMultilinearPCSConfig<F>) -> Self {
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

        let pcs = TensorMultilinearPCS::new(config);

        Self { r1cs, pcs }
    }

    pub fn prove(
        &self,
        r1cs_witness: &[F],
        r1cs_input: &[F],
        transcript: &mut Transcript<F>,
    ) -> (PartialSpartanProof<F>, Vec<F>) {
        // Multilinear extension requires the number of evaluations
        // to be a power of two to uniquely determine the polynomial
        let mut padded_r1cs_witness = r1cs_witness.to_vec();
        padded_r1cs_witness.resize(padded_r1cs_witness.len().next_power_of_two(), F::ZERO);
        let witness_poly = MlPoly::new(padded_r1cs_witness.clone());

        let Z = R1CS::construct_z(r1cs_witness, r1cs_input);

        // Commit the witness polynomial
        let comm_witness_timer = start_timer!(|| "Commit witness");
        let committed_witness = self.pcs.commit(&padded_r1cs_witness);
        let witness_comm = committed_witness.committed_tree.root;
        end_timer!(comm_witness_timer);

        // Add the witness commitment to the transcript
        transcript.append_bytes(&witness_comm);

        // ############################
        // Phase 1
        // ###################

        let m = (self.r1cs.z_len() as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);

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

        let rx = transcript.challenge_vec(m);

        let sc_phase_1_timer = start_timer!(|| "Sumcheck phase 1");

        let sc_phase_1 = SumCheckPhase1::new(
            Az_poly.clone(),
            Bz_poly.clone(),
            Cz_poly.clone(),
            tau.clone(),
            rx.clone(),
        );
        let (sc_proof_1, (v_A, v_B, v_C)) = sc_phase_1.prove(&self.pcs, transcript);
        end_timer!(sc_phase_1_timer);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        // Phase 2
        let r = transcript.challenge_vec(3);

        // T_2 should equal teh evaluations of the random linear combined polynomials

        let ry = transcript.challenge_vec(m);
        let sc_phase_2_timer = start_timer!(|| "Sumcheck phase 2");
        let sc_phase_2 = SumCheckPhase2::new(
            self.r1cs.A.clone(),
            self.r1cs.B.clone(),
            self.r1cs.C.clone(),
            Z.clone(),
            rx.clone(),
            r.as_slice().try_into().unwrap(),
            ry.clone(),
        );

        let sc_proof_2 = sc_phase_2.prove(&self.pcs, transcript);
        end_timer!(sc_phase_2_timer);

        let z_open_timer = start_timer!(|| "Open witness poly");
        // Prove the evaluation of the polynomial Z(y) at ry
        let z_eval_proof = self.pcs.open(
            &committed_witness,
            &padded_r1cs_witness,
            &ry[1..],
            witness_poly.eval(&ry[1..]),
            transcript,
        );
        end_timer!(z_open_timer);

        // Prove the evaluation of the polynomials A(y), B(y), C(y) at ry

        let rx_ry = vec![ry, rx].concat();
        (
            PartialSpartanProof {
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

    pub fn verify_partial(
        &self,
        partial_proof: &PartialSpartanProof<F>,
        transcript: &mut Transcript<F>,
    ) {
        partial_proof.z_comm.append_to_transcript(transcript);

        let A_mle = self.r1cs.A.to_ml_extension();
        let B_mle = self.r1cs.B.to_ml_extension();
        let C_mle = self.r1cs.C.to_ml_extension();

        let m = (self.r1cs.z_len() as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);
        let rx = transcript.challenge_vec(m);

        transcript.append_fe(&partial_proof.sc_proof_1.blinder_poly_sum);
        transcript.append_bytes(&partial_proof.sc_proof_1.blinder_poly_eval_proof.u_hat_comm);

        let rho = transcript.challenge_fe();

        let ex = SumCheckPhase1::verify_round_polys(&partial_proof.sc_proof_1, &rx, rho);

        self.pcs.verify(
            &partial_proof.sc_proof_1.blinder_poly_eval_proof,
            transcript,
        );

        // The final eval should equal
        let v_A = partial_proof.v_A;
        let v_B = partial_proof.v_B;
        let v_C = partial_proof.v_C;

        let T_1_eq = EqPoly::new(tau);
        let T_1 = (v_A * v_B - v_C) * T_1_eq.eval(&rx)
            + rho * partial_proof.sc_proof_1.blinder_poly_eval_proof.y;
        assert_eq!(T_1, ex);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        let r = transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = transcript.challenge_vec(m);

        transcript.append_fe(&partial_proof.sc_proof_2.blinder_poly_sum);
        transcript.append_bytes(&partial_proof.sc_proof_2.blinder_poly_eval_proof.u_hat_comm);

        let rho_2 = transcript.challenge_fe();

        let T_2 =
            (r_A * v_A + r_B * v_B + r_C * v_C) + rho_2 * partial_proof.sc_proof_2.blinder_poly_sum;
        let final_poly_eval =
            SumCheckPhase2::verify_round_polys(T_2, &partial_proof.sc_proof_2, &ry);

        self.pcs.verify(
            &partial_proof.sc_proof_2.blinder_poly_eval_proof,
            transcript,
        );

        assert_eq!(partial_proof.z_eval_proof.x, ry[1..]);
        let rx_ry = [rx, ry.clone()].concat();

        let witness_eval = partial_proof.z_eval_proof.y;
        let A_eval = A_mle.eval(&rx_ry);
        let B_eval = B_mle.eval(&rx_ry);
        let C_eval = C_mle.eval(&rx_ry);

        self.pcs.verify(&partial_proof.z_eval_proof, transcript);

        let input = (0..self.r1cs.num_input)
            .map(|i| (i + 1, partial_proof.pub_input[i]))
            .collect::<Vec<(usize, F)>>();

        let input_poly = SparseMLPoly::new(vec![vec![(0, F::ONE)], input].concat(), ry.len() - 1);
        let input_poly_eval = input_poly.eval(&ry[1..]);

        let z_eval = (F::ONE - ry[0]) * input_poly_eval + ry[0] * witness_eval;

        let T_opened = (r_A * A_eval + r_B * B_eval + r_C * C_eval) * z_eval
            + rho_2 * partial_proof.sc_proof_2.blinder_poly_eval_proof.y;
        assert_eq!(T_opened, final_poly_eval);
    }
}

#[cfg(test)]
mod tests {

    use crate::tensor_pcs::det_num_cols;

    use super::*;

    #[test]
    fn test_shockwave_plus() {
        type F = ark_secp256k1::Fq;

        let num_vars = 20;
        let num_input = 3;
        let l = 2;

        let (r1cs, witness, pub_input) = R1CS::<F>::produce_synthetic_r1cs(num_vars, num_input);

        let config = TensorRSMultilinearPCSConfig::new(r1cs.z_len(), 2, l);

        let ShockwavePlus = ShockwavePlus::new(r1cs.clone(), config);
        let mut prover_transcript = Transcript::new(b"test");

        let (partial_proof, _) = ShockwavePlus.prove(&witness, &pub_input, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"test");
        ShockwavePlus.verify_partial(&partial_proof, &mut verifier_transcript);
    }
}
