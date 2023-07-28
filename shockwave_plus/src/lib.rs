#![allow(non_snake_case)]
mod polynomial;
mod r1cs;
mod sumcheck;
mod utils;

use ark_std::{end_timer, start_timer};
use serde::{Deserialize, Serialize};
use sumcheck::{SCPhase1Proof, SCPhase2Proof, SumCheckPhase1, SumCheckPhase2};

// Exports
pub use r1cs::R1CS;
pub use tensor_pcs::*;

#[derive(Serialize, Deserialize)]
pub struct PartialSpartanProof<F: FieldExt> {
    pub z_comm: [u8; 32],
    pub sc_proof_1: SCPhase1Proof<F>,
    pub sc_proof_2: SCPhase2Proof<F>,
    pub z_eval_proof: TensorMLOpening<F>,
    pub v_A: F,
    pub v_B: F,
    pub v_C: F,
}

pub struct FullSpartanProof<F: FieldExt> {
    pub partial_proof: PartialSpartanProof<F>,
    pub A_eval_proof: TensorMLOpening<F>,
    pub B_eval_proof: TensorMLOpening<F>,
    pub C_eval_proof: TensorMLOpening<F>,
}

pub struct ShockwavePlus<F: FieldExt> {
    pub r1cs: R1CS<F>,
    pub pcs_witness: TensorMultilinearPCS<F>,
}

impl<F: FieldExt> ShockwavePlus<F> {
    pub fn new(r1cs: R1CS<F>, l: usize, num_rows: usize) -> Self {
        let num_cols = r1cs.num_vars / num_rows;

        // Make sure that there are enough columns to run the l queries
        assert!(num_cols > l);

        let expansion_factor = 2;

        let ecfft_config = rs_config::ecfft::gen_config(num_cols);

        let pcs_config = TensorRSMultilinearPCSConfig::<F> {
            expansion_factor,
            domain_powers: None,
            fft_domain: None,
            ecfft_config: Some(ecfft_config),
            l,
            num_entries: r1cs.num_vars,
            num_rows,
        };

        let pcs_witness = TensorMultilinearPCS::new(pcs_config);
        Self { r1cs, pcs_witness }
    }

    pub fn prove(
        &self,
        witness: &[F],
        transcript: &mut Transcript<F>,
    ) -> (PartialSpartanProof<F>, Vec<F>) {
        // Compute the multilinear extension of the witness
        assert!(witness.len().is_power_of_two());
        let witness_poly = SparseMLPoly::from_dense(witness.to_vec());

        // Commit the witness polynomial
        let comm_witness_timer = start_timer!(|| "Commit witness");
        let committed_witness = self.pcs_witness.commit(&witness_poly);
        let witness_comm = committed_witness.committed_tree.root;
        end_timer!(comm_witness_timer);

        transcript.append_bytes(&witness_comm);

        // ############################
        // Phase 1: The sum-checks
        // ###################

        let m = (self.r1cs.num_vars as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);
        let mut tau_rev = tau.clone();
        tau_rev.reverse();

        // First
        // Compute the multilinear extension of the R1CS matrices.
        // Prove that he Q_poly is a zero-polynomial

        // Q_poly is a zero-polynomial iff F_io evaluates to zero
        // over the m-dimensional boolean hypercube..

        // We prove using the sum-check protocol.

        // G_poly = A_poly * B_poly - C_poly

        let num_rows = self.r1cs.num_cons;
        let Az_poly = self.r1cs.A.mul_vector(num_rows, witness);
        let Bz_poly = self.r1cs.B.mul_vector(num_rows, witness);
        let Cz_poly = self.r1cs.C.mul_vector(num_rows, witness);

        // Prove that the polynomial Q(t)
        // \sum_{x \in {0, 1}^m} (Az_poly(x) * Bz_poly(x) - Cz_poly(x)) eq(tau, x)
        // is a zero-polynomial using the sum-check protocol.

        let rx = transcript.challenge_vec(m);
        let mut rx_rev = rx.clone();
        rx_rev.reverse();

        let sc_phase_1_timer = start_timer!(|| "Sumcheck phase 1");

        let sc_phase_1 = SumCheckPhase1::new(
            Az_poly.clone(),
            Bz_poly.clone(),
            Cz_poly.clone(),
            tau_rev.clone(),
            rx.clone(),
        );
        let (sc_proof_1, (v_A, v_B, v_C)) = sc_phase_1.prove(transcript);
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
            witness.to_vec(),
            rx.clone(),
            r.as_slice().try_into().unwrap(),
            ry.clone(),
        );

        let sc_proof_2 = sc_phase_2.prove(transcript);
        end_timer!(sc_phase_2_timer);

        let mut ry_rev = ry.clone();
        ry_rev.reverse();
        let z_open_timer = start_timer!(|| "Open witness poly");
        // Prove the evaluation of the polynomial Z(y) at ry
        let z_eval_proof =
            self.pcs_witness
                .open(&committed_witness, &witness_poly, &ry_rev, transcript);
        end_timer!(z_open_timer);

        // Prove the evaluation of the polynomials A(y), B(y), C(y) at ry

        let rx_ry = vec![ry_rev, rx_rev].concat();
        (
            PartialSpartanProof {
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

        let m = (self.r1cs.num_vars as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);
        let rx = transcript.challenge_vec(m);
        let mut rx_rev = rx.clone();
        rx_rev.reverse();

        transcript.append_fe(&partial_proof.sc_proof_1.blinder_poly_sum);
        let rho = transcript.challenge_fe();

        let ex = SumCheckPhase1::verify_round_polys(&partial_proof.sc_proof_1, &rx, rho);

        // The final eval should equal
        let v_A = partial_proof.v_A;
        let v_B = partial_proof.v_B;
        let v_C = partial_proof.v_C;

        let T_1_eq = EqPoly::new(tau);
        let T_1 = (v_A * v_B - v_C) * T_1_eq.eval(&rx_rev)
            + rho * partial_proof.sc_proof_1.blinder_poly_eval_claim;
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
        let rho_2 = transcript.challenge_fe();

        let T_2 =
            (r_A * v_A + r_B * v_B + r_C * v_C) + rho_2 * partial_proof.sc_proof_2.blinder_poly_sum;
        let final_poly_eval =
            SumCheckPhase2::verify_round_polys(T_2, &partial_proof.sc_proof_2, &ry);

        let mut ry_rev = ry.clone();
        ry_rev.reverse();

        let rx_ry = [rx, ry].concat();
        assert_eq!(partial_proof.z_eval_proof.x, ry_rev);

        let z_eval = partial_proof.z_eval_proof.y;
        let A_eval = A_mle.eval(&rx_ry);
        let B_eval = B_mle.eval(&rx_ry);
        let C_eval = C_mle.eval(&rx_ry);

        self.pcs_witness.verify(
            &partial_proof.z_eval_proof,
            &partial_proof.z_comm,
            transcript,
        );

        let T_opened = (r_A * A_eval + r_B * B_eval + r_C * C_eval) * z_eval
            + rho_2 * partial_proof.sc_proof_2.blinder_poly_eval_claim;
        assert_eq!(T_opened, final_poly_eval);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shockwave_plus() {
        type F = halo2curves::secp256k1::Fp;

        let num_cons = 2usize.pow(6);
        let num_vars = num_cons;
        let num_input = 0;
        let l = 10;

        let (r1cs, witness) = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        let num_rows = 4;
        let ShockwavePlus = ShockwavePlus::new(r1cs.clone(), l, num_rows);
        let mut prover_transcript = Transcript::new(b"bench");
        let (partial_proof, _) = ShockwavePlus.prove(&witness, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"bench");
        ShockwavePlus.verify_partial(&partial_proof, &mut verifier_transcript);
    }
}
