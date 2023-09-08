use crate::polynomial::{eq_poly::EqPoly, ml_poly::MlPoly};
use crate::sumcheck::unipoly::UniPoly;
use crate::sumcheck::SumCheckProof;
use crate::tensor_pcs::TensorMultilinearPCS;
use crate::transcript::Transcript;

use crate::FieldGC;

pub struct SumCheckPhase1<F: FieldGC> {
    Az_evals: Vec<F>,
    Bz_evals: Vec<F>,
    Cz_evals: Vec<F>,
    bound_eq_poly: EqPoly<F>,
    challenge: Vec<F>,
}

impl<F: FieldGC> SumCheckPhase1<F> {
    pub fn new(
        Az_evals: Vec<F>,
        Bz_evals: Vec<F>,
        Cz_evals: Vec<F>,
        tau: Vec<F>,
        challenge: Vec<F>,
    ) -> Self {
        let bound_eq_poly = EqPoly::new(tau);
        Self {
            Az_evals,
            Bz_evals,
            Cz_evals,
            bound_eq_poly,
            challenge,
        }
    }

    // TODO: DRY prove and prove_zk

    pub fn prove(&self) -> (SumCheckProof<F>, (F, F, F)) {
        let num_vars = (self.Az_evals.len() as f64).log2() as usize;
        let mut round_poly_coeffs = Vec::<Vec<F>>::with_capacity(num_vars - 1);

        let mut A_table = self.Az_evals.clone();
        let mut B_table = self.Bz_evals.clone();
        let mut C_table = self.Cz_evals.clone();
        let mut eq_table = self.bound_eq_poly.evals();

        let zero = F::ZERO;
        let one = F::ONE;
        let two = F::from(2u64);
        let three = F::from(3u64);

        for j in 0..num_vars {
            let r_i = self.challenge[j];

            let high_index = 2usize.pow((num_vars - j - 1) as u32);

            let mut evals = [F::ZERO; 4];

            // https://eprint.iacr.org/2019/317.pdf#subsection.3.2
            for b in 0..high_index {
                for (i, eval_at) in [zero, one, two, three].iter().enumerate() {
                    let a_eval = A_table[b] + (A_table[b + high_index] - A_table[b]) * eval_at;
                    let b_eval = B_table[b] + (B_table[b + high_index] - B_table[b]) * eval_at;
                    let c_eval = C_table[b] + (C_table[b + high_index] - C_table[b]) * eval_at;
                    let eq_eval = eq_table[b] + (eq_table[b + high_index] - eq_table[b]) * eval_at;
                    evals[i] += (a_eval * b_eval - c_eval) * eq_eval;
                }

                A_table[b] = A_table[b] + (A_table[b + high_index] - A_table[b]) * r_i;
                B_table[b] = B_table[b] + (B_table[b + high_index] - B_table[b]) * r_i;
                C_table[b] = C_table[b] + (C_table[b + high_index] - C_table[b]) * r_i;
                eq_table[b] = eq_table[b] + (eq_table[b + high_index] - eq_table[b]) * r_i;
            }

            let round_poly = UniPoly::interpolate(&evals);

            round_poly_coeffs.push(round_poly.coeffs);
        }

        let v_A = A_table[0];
        let v_B = B_table[0];
        let v_C = C_table[0];

        (
            SumCheckProof {
                round_poly_coeffs,
                blinder_poly_eval_proof: None,
                blinder_poly_sum: None,
            },
            (v_A, v_B, v_C),
        )
    }

    pub fn prove_zk(
        &self,
        pcs: &TensorMultilinearPCS<F>,
        transcript: &mut Transcript<F>,
    ) -> (SumCheckProof<F>, (F, F, F)) {
        let num_vars = (self.Az_evals.len() as f64).log2() as usize;
        let mut round_poly_coeffs = Vec::<Vec<F>>::with_capacity(num_vars - 1);

        // We implement the zero-knowledge sumcheck protocol
        // described in Section 4.1 https://eprint.iacr.org/2019/317.pdf

        let mut rng = rand::thread_rng();
        // Sample a blinding polynomial g(x_1, ..., x_m)
        let blinder_poly_evals = (0..2usize.pow(num_vars as u32))
            .map(|_| F::rand(&mut rng))
            .collect::<Vec<F>>();
        let blinder_poly = MlPoly::new(blinder_poly_evals.clone());
        let blinder_poly_sum = blinder_poly_evals.iter().fold(F::ZERO, |acc, x| acc + x);

        let blinder_poly_comm = pcs.commit(&blinder_poly_evals, true);

        transcript.append_fe(&blinder_poly_sum);
        transcript.append_bytes(&blinder_poly_comm.committed_tree.root);

        let rho = transcript.challenge_fe();

        // Compute the sum of g(x_1, ... x_m) over the boolean hypercube

        // Do the sum check for f + \rho g

        let mut A_table = self.Az_evals.clone();
        let mut B_table = self.Bz_evals.clone();
        let mut C_table = self.Cz_evals.clone();
        let mut blinder_table = blinder_poly_evals.clone();
        let mut eq_table = self.bound_eq_poly.evals();

        let zero = F::ZERO;
        let one = F::ONE;
        let two = F::from(2u64);
        let three = F::from(3u64);

        for j in 0..num_vars {
            let r_i = self.challenge[j];

            let high_index = 2usize.pow((num_vars - j - 1) as u32);

            let mut evals = [F::ZERO; 4];

            // https://eprint.iacr.org/2019/317.pdf#subsection.3.2
            for b in 0..high_index {
                for (i, eval_at) in [zero, one, two, three].iter().enumerate() {
                    let a_eval = A_table[b] + (A_table[b + high_index] - A_table[b]) * eval_at;
                    let b_eval = B_table[b] + (B_table[b + high_index] - B_table[b]) * eval_at;
                    let c_eval = C_table[b] + (C_table[b + high_index] - C_table[b]) * eval_at;
                    let eq_eval = eq_table[b] + (eq_table[b + high_index] - eq_table[b]) * eval_at;
                    let blinder_eval = blinder_table[b]
                        + (blinder_table[b + high_index] - blinder_table[b]) * eval_at;
                    evals[i] += ((a_eval * b_eval - c_eval) * eq_eval) + rho * blinder_eval;
                }

                A_table[b] = A_table[b] + (A_table[b + high_index] - A_table[b]) * r_i;
                B_table[b] = B_table[b] + (B_table[b + high_index] - B_table[b]) * r_i;
                C_table[b] = C_table[b] + (C_table[b + high_index] - C_table[b]) * r_i;
                eq_table[b] = eq_table[b] + (eq_table[b + high_index] - eq_table[b]) * r_i;
                blinder_table[b] =
                    blinder_table[b] + (blinder_table[b + high_index] - blinder_table[b]) * r_i;
            }

            // TODO: Maybe send the evaluations to the verifier?
            let round_poly = UniPoly::interpolate(&evals);

            round_poly_coeffs.push(round_poly.coeffs);
        }

        let v_A = A_table[0];
        let v_B = B_table[0];
        let v_C = C_table[0];

        // Prove the evaluation of the blinder polynomial at rx.
        let blinder_poly_eval_proof = pcs.open(
            &blinder_poly_comm,
            &blinder_poly_evals,
            &self.challenge,
            blinder_poly.eval(&self.challenge),
            transcript,
            true,
        );

        (
            SumCheckProof {
                round_poly_coeffs,
                blinder_poly_eval_proof: Some(blinder_poly_eval_proof),
                blinder_poly_sum: Some(blinder_poly_sum),
            },
            (v_A, v_B, v_C),
        )
    }

    pub fn verify_round_polys(proof: &SumCheckProof<F>, challenge: &[F], rho: Option<F>) -> F {
        debug_assert_eq!(proof.round_poly_coeffs.len(), challenge.len());

        let mut target = if proof.is_blinded() {
            rho.unwrap() * proof.blinder_poly_sum.unwrap()
        } else {
            F::ZERO
        };

        for (i, coeffs) in proof.round_poly_coeffs.iter().enumerate() {
            let round_poly = UniPoly::new(coeffs.clone());
            assert_eq!(
                round_poly.eval(F::ZERO) + round_poly.eval(F::ONE),
                target,
                "round poly {} failed",
                i
            );
            target = round_poly.eval(challenge[i]);
        }

        target
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = ark_secp256k1::Fq;
    use ark_ff::Field;

    #[test]
    fn test_unipoly_3() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64), F::from(4u64)];
        let eval_at = F::from(33);

        let mut expected_eval = F::ZERO;
        for i in 0..coeffs.len() {
            expected_eval += coeffs[i] * eval_at.pow(&[3 - i as u64, 0, 0, 0]);
        }

        let mut evals = [F::ZERO; 4];
        for i in 0..4 {
            let eval_at = F::from(i as u64);
            let mut eval_i = F::ZERO;
            for j in 0..coeffs.len() {
                eval_i += coeffs[j] * eval_at.pow(&[3 - j as u64, 0, 0, 0]);
            }
            evals[i] = eval_i;
        }

        let uni_poly = UniPoly::interpolate(&evals);
        let eval = uni_poly.eval(eval_at);
        assert_eq!(eval, expected_eval);
    }
}
