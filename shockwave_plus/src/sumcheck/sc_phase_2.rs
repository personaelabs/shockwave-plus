use crate::r1cs::r1cs::Matrix;
use crate::sumcheck::unipoly::UniPoly;
use crate::FieldExt;
use serde::{Deserialize, Serialize};
use tensor_pcs::{EqPoly, SparseMLPoly, TensorMLOpening, TensorMultilinearPCS, Transcript};

#[derive(Serialize, Deserialize)]
pub struct SCPhase2Proof<F: FieldExt> {
    pub round_polys: Vec<UniPoly<F>>,
    pub blinder_poly_sum: F,
    pub blinder_poly_eval_proof: TensorMLOpening<F>,
}

pub struct SumCheckPhase2<F: FieldExt> {
    A_mat: Matrix<F>,
    B_mat: Matrix<F>,
    C_mat: Matrix<F>,
    Z_evals: Vec<F>,
    rx: Vec<F>,
    r: [F; 3],
    challenge: Vec<F>,
}

impl<F: FieldExt> SumCheckPhase2<F> {
    pub fn new(
        A_mat: Matrix<F>,
        B_mat: Matrix<F>,
        C_mat: Matrix<F>,
        Z_evals: Vec<F>,
        rx: Vec<F>,
        r: [F; 3],
        challenge: Vec<F>,
    ) -> Self {
        Self {
            A_mat,
            B_mat,
            C_mat,
            Z_evals,
            rx,
            r,
            challenge,
        }
    }

    pub fn prove(
        &self,
        pcs: &TensorMultilinearPCS<F>,
        transcript: &mut Transcript<F>,
    ) -> SCPhase2Proof<F> {
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

        let mut rng = rand::thread_rng();
        // Sample a blinding polynomial g(x_1, ..., x_m) of degree 3
        let random_evals = (0..2usize.pow(num_vars as u32))
            .map(|_| F::random(&mut rng))
            .collect::<Vec<F>>();
        let blinder_poly_sum = random_evals.iter().fold(F::ZERO, |acc, x| acc + x);
        let blinder_poly = SparseMLPoly::from_dense(random_evals);
        let blinder_poly_comm = pcs.commit(&blinder_poly);

        transcript.append_fe(&blinder_poly_sum);
        transcript.append_bytes(&blinder_poly_comm.committed_tree.root);

        let rho = transcript.challenge_fe();

        let mut round_polys: Vec<UniPoly<F>> = Vec::<UniPoly<F>>::with_capacity(num_vars);

        let mut A_table = A_evals.clone();
        let mut B_table = B_evals.clone();
        let mut C_table = C_evals.clone();
        let mut Z_table = self.Z_evals.clone();
        let mut blinder_table = blinder_poly
            .evals
            .iter()
            .map(|(_, x)| *x)
            .collect::<Vec<F>>();

        let zero = F::ZERO;
        let one = F::ONE;
        let two = F::from(2);

        for j in 0..num_vars {
            let high_index = 2usize.pow((num_vars - j - 1) as u32);
            let mut evals = [F::ZERO; 3];

            for b in 0..high_index {
                let r_y_i = self.challenge[j];
                for (i, eval_at) in [zero, one, two].iter().enumerate() {
                    let a_eval = A_table[b] + (A_table[b + high_index] - A_table[b]) * eval_at;
                    let b_eval = B_table[b] + (B_table[b + high_index] - B_table[b]) * eval_at;
                    let c_eval = C_table[b] + (C_table[b + high_index] - C_table[b]) * eval_at;
                    let z_eval = Z_table[b] + (Z_table[b + high_index] - Z_table[b]) * eval_at;
                    let blinder_eval = blinder_table[b]
                        + (blinder_table[b + high_index] - blinder_table[b]) * eval_at;
                    evals[i] +=
                        (a_eval * r_A + b_eval * r_B + c_eval * r_C) * z_eval + rho * blinder_eval;
                }

                A_table[b] = A_table[b] + (A_table[b + high_index] - A_table[b]) * r_y_i;
                B_table[b] = B_table[b] + (B_table[b + high_index] - B_table[b]) * r_y_i;
                C_table[b] = C_table[b] + (C_table[b + high_index] - C_table[b]) * r_y_i;
                Z_table[b] = Z_table[b] + (Z_table[b + high_index] - Z_table[b]) * r_y_i;
                blinder_table[b] =
                    blinder_table[b] + (blinder_table[b + high_index] - blinder_table[b]) * r_y_i;
            }

            let round_poly = UniPoly::interpolate(&evals);
            round_polys.push(round_poly);
        }

        let ry = self.challenge.clone();

        let blinder_poly_eval_proof = pcs.open(&blinder_poly_comm, &blinder_poly, &ry, transcript);

        SCPhase2Proof {
            round_polys,
            blinder_poly_sum,
            blinder_poly_eval_proof,
        }
    }

    pub fn verify_round_polys(sum_target: F, proof: &SCPhase2Proof<F>, challenge: &[F]) -> F {
        debug_assert_eq!(proof.round_polys.len(), challenge.len());

        let zero = F::ZERO;
        let one = F::ONE;

        let mut target = sum_target;
        for (i, round_poly) in proof.round_polys.iter().enumerate() {
            assert_eq!(
                round_poly.eval(zero) + round_poly.eval(one),
                target,
                "i = {}",
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
    use halo2curves::ff::Field;
    use halo2curves::secp256k1::Fp;
    type F = Fp;

    #[test]
    fn test_unipoly_2() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64)];
        let eval_at = Fp::from(33);

        let mut expected_eval = F::ZERO;
        for i in 0..coeffs.len() {
            expected_eval += coeffs[i] * eval_at.pow(&[i as u64, 0, 0, 0]);
        }

        let mut evals = [F::ZERO; 3];
        for i in 0..3 {
            let eval_at = F::from(i as u64);
            let mut eval_i = F::ZERO;
            for j in 0..coeffs.len() {
                eval_i += coeffs[j] * eval_at.pow(&[j as u64, 0, 0, 0]);
            }
            evals[i] = eval_i;
        }

        let uni_poly = UniPoly::interpolate(&evals);
        let eval = uni_poly.eval(eval_at);
        assert_eq!(eval, expected_eval);
    }
}
