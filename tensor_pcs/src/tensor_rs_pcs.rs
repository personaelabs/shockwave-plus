use crate::rs_config::ecfft::ECFFTConfig;
use crate::tree::BaseOpening;
use crate::FieldExt;
use ecfft::extend;
use serde::{Deserialize, Serialize};

use crate::polynomial::eq_poly::EqPoly;
use crate::tensor_code::TensorCode;
use crate::transcript::Transcript;
use crate::utils::{det_num_cols, det_num_rows, dot_prod, hash_all, rlc_rows, sample_indices};

use super::tensor_code::CommittedTensorCode;

#[derive(Clone)]
pub struct TensorRSMultilinearPCSConfig<F: FieldExt> {
    pub expansion_factor: usize,
    pub ecfft_config: ECFFTConfig<F>,
    pub l: usize,
}

impl<F: FieldExt> TensorRSMultilinearPCSConfig<F> {
    pub fn num_cols(&self, num_entries: usize) -> usize {
        det_num_cols(num_entries, self.l)
    }

    pub fn num_rows(&self, num_entries: usize) -> usize {
        det_num_rows(num_entries, self.l)
    }
}

#[derive(Clone)]
pub struct TensorMultilinearPCS<F: FieldExt> {
    config: TensorRSMultilinearPCSConfig<F>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct TensorMLOpening<F: FieldExt> {
    pub x: Vec<F>,
    pub y: F,
    pub base_opening: BaseOpening,
    pub test_query_leaves: Vec<Vec<F>>,
    pub eval_query_leaves: Vec<Vec<F>>,
    pub u_hat_comm: [u8; 32],
    pub test_u_prime: Vec<F>,
    pub test_r_prime: Vec<F>,
    pub eval_r_prime: Vec<F>,
    pub eval_u_prime: Vec<F>,
    pub poly_num_vars: usize,
}

impl<F: FieldExt> TensorMultilinearPCS<F> {
    pub fn new(config: TensorRSMultilinearPCSConfig<F>) -> Self {
        Self { config }
    }

    pub fn commit(&self, ml_poly_evals: &[F]) -> CommittedTensorCode<F> {
        // Merkle commit to the evaluations of the polynomial
        let n = ml_poly_evals.len();
        assert!(n.is_power_of_two());
        let tensor_code = self.encode_zk(ml_poly_evals);
        let tree = tensor_code.commit(self.config.num_cols(n), self.config.num_rows(n));
        tree
    }

    pub fn open(
        &self,
        u_hat_comm: &CommittedTensorCode<F>,
        // TODO: Remove poly and use u_hat_comm
        ml_poly_evals: &[F],
        point: &[F],
        eval: F,
        transcript: &mut Transcript<F>,
    ) -> TensorMLOpening<F> {
        let n = ml_poly_evals.len();
        assert!(n.is_power_of_two());
        let num_vars = (n as f64).log2() as usize;

        let num_cols = self.config.num_cols(n);
        let num_rows = self.config.num_rows(n);

        debug_assert_eq!(num_vars, point.len());

        // ########################################
        // Testing phase
        // Prove the consistency between the random linear combination of the evaluation tensor (u_prime)
        // and the tensor codeword (u_hat)
        // ########################################

        // Derive the challenge vector;
        let r_u = transcript.challenge_vec(num_rows);

        let u = (0..num_rows)
            .map(|i| ml_poly_evals[(i * num_cols)..((i + 1) * num_cols)].to_vec())
            .collect::<Vec<Vec<F>>>();

        // Random linear combination of the rows of the polynomial in a tensor structure
        let test_u_prime = rlc_rows(u.clone(), &r_u);

        // Random linear combination of the blinder
        let blinder = u_hat_comm
            .tensor_codeword
            .0
            .iter()
            .map(|row| row[(row.len() / 2)..].to_vec())
            .collect::<Vec<Vec<F>>>();

        debug_assert_eq!(blinder[0].len(), u_hat_comm.tensor_codeword.0[0].len() / 2);

        let test_r_prime = rlc_rows(blinder.clone(), &r_u);

        let num_indices = self.config.l;
        let indices = sample_indices(num_indices, num_cols * 2, transcript);

        let test_queries = self.test_phase(&indices, &u_hat_comm);

        // ########################################
        // Evaluation phase
        // Prove the consistency
        // ########################################

        // Get the evaluation point

        let log2_num_rows = (num_rows as f64).log2() as usize;
        let q1 = EqPoly::new(point[0..log2_num_rows].to_vec()).evals();

        let eval_r_prime = rlc_rows(blinder, &q1);

        let eval_u_prime = rlc_rows(u.clone(), &q1);

        let eval_queries = self.test_phase(&indices, &u_hat_comm);

        TensorMLOpening {
            x: point.to_vec(),
            y: eval,
            eval_query_leaves: eval_queries,
            test_query_leaves: test_queries,
            u_hat_comm: u_hat_comm.committed_tree.root(),
            test_u_prime,
            test_r_prime,
            eval_r_prime,
            eval_u_prime,
            base_opening: BaseOpening {
                hashes: u_hat_comm.committed_tree.column_roots.clone(),
            },
            poly_num_vars: num_vars,
        }
    }
}

impl<F: FieldExt> TensorMultilinearPCS<F> {
    pub fn verify(&self, opening: &TensorMLOpening<F>, transcript: &mut Transcript<F>) {
        let poly_num_entries = 2usize.pow(opening.poly_num_vars as u32);
        let num_rows = self.config.num_rows(poly_num_entries);
        let num_cols = self.config.num_cols(poly_num_entries);

        // Verify the base opening
        let base_opening = &opening.base_opening;
        base_opening.verify(opening.u_hat_comm);

        // ########################################
        // Verify test phase
        // ########################################

        let r_u = transcript.challenge_vec(num_rows);

        let test_u_prime_rs_codeword = self
            .rs_encode(&opening.test_u_prime)
            .iter()
            .zip(opening.test_r_prime.iter())
            .map(|(c, r)| *c + *r)
            .collect::<Vec<F>>();

        let num_indices = self.config.l;
        let indices = sample_indices(num_indices, num_cols * 2, transcript);

        debug_assert_eq!(indices.len(), opening.test_query_leaves.len());
        for (expected_index, leaves) in indices.iter().zip(opening.test_query_leaves.iter()) {
            // Verify that the hashes of the leaves equals the corresponding column root
            let leaf_bytes = leaves
                .iter()
                .map(|x| x.to_repr())
                .collect::<Vec<[u8; 32]>>();
            let column_root = hash_all(&leaf_bytes);
            let expected_column_root = base_opening.hashes[*expected_index];
            assert_eq!(column_root, expected_column_root);

            let mut sum = F::ZERO;
            for (leaf, r_i) in leaves.iter().zip(r_u.iter()) {
                sum += *r_i * *leaf;
            }
            assert_eq!(sum, test_u_prime_rs_codeword[*expected_index]);
        }

        // ########################################
        // Verify evaluation phase
        // ########################################

        let log2_num_rows = (num_rows as f64).log2() as usize;
        let q1 = EqPoly::new(opening.x[0..log2_num_rows].to_vec()).evals();
        let q2 = EqPoly::new(opening.x[log2_num_rows..].to_vec()).evals();

        let eval_u_prime_rs_codeword = self
            .rs_encode(&opening.eval_u_prime)
            .iter()
            .zip(opening.eval_r_prime.iter())
            .map(|(c, r)| *c + *r)
            .collect::<Vec<F>>();

        debug_assert_eq!(q1.len(), opening.eval_query_leaves[0].len());
        debug_assert_eq!(indices.len(), opening.test_query_leaves.len());
        for (expected_index, leaves) in indices.iter().zip(opening.eval_query_leaves.iter()) {
            // TODO: Don't need to check the leaves again?
            // Verify that the hashes of the leaves equals the corresponding column root
            let leaf_bytes = leaves
                .iter()
                .map(|x| x.to_repr())
                .collect::<Vec<[u8; 32]>>();
            let column_root = hash_all(&leaf_bytes);
            let expected_column_root = base_opening.hashes[*expected_index];
            assert_eq!(column_root, expected_column_root);

            let mut sum = F::ZERO;
            for (leaf, q1_i) in leaves.iter().zip(q1.iter()) {
                sum += *q1_i * *leaf;
            }
            assert_eq!(sum, eval_u_prime_rs_codeword[*expected_index]);
        }

        let expected_eval = dot_prod(&opening.eval_u_prime, &q2);
        assert_eq!(expected_eval, opening.y);
    }

    fn split_encode(&self, message: &[F]) -> Vec<F> {
        let codeword = self.rs_encode(message);

        let mut rng = rand::thread_rng();
        let blinder = (0..codeword.len())
            .map(|_| F::random(&mut rng))
            .collect::<Vec<F>>();

        let mut randomized_codeword = codeword
            .iter()
            .zip(blinder.clone().iter())
            .map(|(c, b)| *b + *c)
            .collect::<Vec<F>>();

        randomized_codeword.extend_from_slice(&blinder);
        debug_assert_eq!(randomized_codeword.len(), codeword.len() * 2);
        randomized_codeword
    }

    fn rs_encode(&self, message: &[F]) -> Vec<F> {
        let mut padded_message = message.to_vec();
        padded_message.resize(message.len().next_power_of_two(), F::ZERO);
        let codeword_len = padded_message.len() * self.config.expansion_factor;

        let codeword_len_log2 = (codeword_len as f64).log2() as usize;

        // Resize the domain to the correct size

        let ecfft_config = &self.config.ecfft_config;
        let config_domain_size = ecfft_config.domain.len();

        assert!(config_domain_size >= codeword_len_log2 - 1);
        let domain = ecfft_config.domain[(config_domain_size - (codeword_len_log2 - 1))..].to_vec();
        let matrices =
            ecfft_config.matrices[(config_domain_size - (codeword_len_log2 - 1))..].to_vec();
        let inverse_matrices = ecfft_config.inverse_matrices
            [(config_domain_size - (codeword_len_log2 - 1))..]
            .to_vec();

        assert_eq!(
            padded_message.len() * self.config.expansion_factor,
            domain[0].len()
        );

        let extended_evals = extend(&padded_message, &domain, &matrices, &inverse_matrices, 0);

        let codeword = [message.to_vec(), extended_evals].concat();

        codeword
    }

    fn test_phase(&self, indices: &[usize], u_hat_comm: &CommittedTensorCode<F>) -> Vec<Vec<F>> {
        // Query the columns of u_hat
        let num_indices = self.config.l;

        let u_hat_openings = indices
            .iter()
            .map(|index| u_hat_comm.query_column(*index))
            .collect::<Vec<Vec<F>>>();

        debug_assert_eq!(u_hat_openings.len(), num_indices);

        u_hat_openings
    }

    fn encode_zk(&self, ml_poly_evals: &[F]) -> TensorCode<F> {
        let n = ml_poly_evals.len();
        assert!(n.is_power_of_two());

        let num_rows = self.config.num_rows(n);
        let num_cols = self.config.num_cols(n);
        debug_assert_eq!(n, num_cols * num_rows);

        let codewords = (0..num_rows)
            .map(|i| &ml_poly_evals[i * num_cols..(i + 1) * num_cols])
            .map(|row| self.split_encode(&row))
            .collect::<Vec<Vec<F>>>();

        TensorCode(codewords)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::ml_poly::MlPoly;
    use crate::rs_config::{ecfft, good_curves::secp256k1::secp256k1_good_curve};

    const TEST_NUM_VARS: usize = 8;
    const TEST_L: usize = 10;

    fn test_poly_evals<F: FieldExt>() -> MlPoly<F> {
        let num_entries: usize = 2usize.pow(TEST_NUM_VARS as u32);

        let evals = (0..num_entries)
            .map(|i| F::from((i + 1) as u64))
            .collect::<Vec<F>>();

        MlPoly::new(evals)
    }

    fn prove_and_verify<F: FieldExt>(ml_poly: &MlPoly<F>, pcs: TensorMultilinearPCS<F>) {
        let ml_poly_evals = &ml_poly.evals;

        let comm = pcs.commit(ml_poly_evals);

        let ml_poly_num_vars = (ml_poly_evals.len() as f64).log2() as usize;
        let open_at = (0..ml_poly_num_vars)
            .map(|i| F::from(i as u64))
            .collect::<Vec<F>>();
        let y = ml_poly.eval(&open_at);

        let mut prover_transcript = Transcript::<F>::new(b"test");
        prover_transcript.append_bytes(&comm.committed_tree.root);
        let opening = pcs.open(&comm, ml_poly_evals, &open_at, y, &mut prover_transcript);

        let mut verifier_transcript = Transcript::<F>::new(b"test");
        verifier_transcript.append_bytes(&comm.committed_tree.root);
        pcs.verify(&opening, &mut verifier_transcript);
    }

    #[test]
    fn test_tensor_pcs() {
        type F = halo2curves::secp256k1::Fp;
        let ml_poly = test_poly_evals();

        let expansion_factor = 2;

        let n = ml_poly.evals.len();
        let num_cols = det_num_cols(n, TEST_L);
        let k = ((num_cols * expansion_factor).next_power_of_two() as f64).log2() as usize;
        let (curve, coset_offset) = secp256k1_good_curve(k);
        let ecfft_config = ecfft::gen_config_form_curve(curve, coset_offset);

        let config = TensorRSMultilinearPCSConfig::<F> {
            expansion_factor,
            ecfft_config,
            l: TEST_L,
        };

        let tensor_pcs_ecf = TensorMultilinearPCS::<F>::new(config);
        prove_and_verify(&ml_poly, tensor_pcs_ecf);
    }
}
