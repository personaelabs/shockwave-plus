use crate::polynomial::ml_poly::MlPoly;
use crate::sumcheck::unipoly::UniPoly;
use serde::{Deserialize, Serialize};
use tensor_pcs::{EqPoly, Transcript};

use crate::FieldExt;

#[derive(Serialize, Deserialize)]
pub struct SCPhase1Proof<F: FieldExt> {
    pub blinder_poly_sum: F,
    pub blinder_poly_eval_claim: F,
    pub round_polys: Vec<UniPoly<F>>,
}

pub struct SumCheckPhase1<F: FieldExt> {
    Az_evals: Vec<F>,
    Bz_evals: Vec<F>,
    Cz_evals: Vec<F>,
    bound_eq_poly: EqPoly<F>,
    challenge: Vec<F>,
}

impl<F: FieldExt> SumCheckPhase1<F> {
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

    pub fn prove(&self, transcript: &mut Transcript<F>) -> (SCPhase1Proof<F>, (F, F, F)) {
        let num_vars = (self.Az_evals.len() as f64).log2() as usize;
        let mut round_polys = Vec::<UniPoly<F>>::with_capacity(num_vars - 1);

        let mut rng = rand::thread_rng();
        // Sample a blinding polynomial g(x_1, ..., x_m) of degree 3
        let random_evals = (0..2usize.pow(num_vars as u32))
            .map(|_| F::random(&mut rng))
            .collect::<Vec<F>>();
        let blinder_poly_sum = random_evals.iter().fold(F::ZERO, |acc, x| acc + x);
        let blinder_poly = MlPoly::new(random_evals);

        transcript.append_fe(&blinder_poly_sum);
        let rho = transcript.challenge_fe();

        // Compute the sum of g(x_1, ... x_m) over the boolean hypercube

        // Do the sum check for f + \rho g

        let mut A_table = self.Az_evals.clone();
        let mut B_table = self.Bz_evals.clone();
        let mut C_table = self.Cz_evals.clone();
        let mut blinder_table = blinder_poly.evals.clone();
        let mut eq_table = self.bound_eq_poly.evals();

        let zero = F::ZERO;
        let one = F::ONE;
        let two = F::from(2);
        let three = F::from(3);

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

            let round_poly = UniPoly::interpolate(&evals);

            round_polys.push(round_poly);
        }

        let v_A = A_table[0];
        let v_B = B_table[0];
        let v_C = C_table[0];

        let rx = self.challenge.clone();
        let blinder_poly_eval_claim = blinder_poly.eval(&rx);

        // Prove the evaluation of the blinder polynomial at rx.

        (
            SCPhase1Proof {
                blinder_poly_sum,
                round_polys,
                blinder_poly_eval_claim,
            },
            (v_A, v_B, v_C),
        )
    }

    pub fn verify_round_polys(proof: &SCPhase1Proof<F>, challenge: &[F], rho: F) -> F {
        debug_assert_eq!(proof.round_polys.len(), challenge.len());

        let zero = F::ZERO;
        let one = F::ONE;

        // target = 0 + rho * blinder_poly_sum
        let mut target = rho * proof.blinder_poly_sum;
        for (i, round_poly) in proof.round_polys.iter().enumerate() {
            assert_eq!(
                round_poly.eval(zero) + round_poly.eval(one),
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
    use halo2curves::secp256k1::Fp;
    type F = Fp;
    use halo2curves::ff::Field;

    #[test]
    fn test_unipoly_3() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64), F::from(4u64)];
        let eval_at = Fp::from(33);

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
