use crate::rs_config::ecfft::ECFFTConfig;
use crate::tree::BaseOpening;
use crate::FieldExt;
use ecfft::extend;
use serde::{Deserialize, Serialize};

use crate::fft::fft;
use crate::polynomial::eq_poly::EqPoly;
use crate::polynomial::sparse_ml_poly::SparseMLPoly;
use crate::tensor_code::TensorCode;
use crate::transcript::Transcript;
use crate::utils::{dot_prod, hash_all, rlc_rows, sample_indices};

use super::tensor_code::CommittedTensorCode;

#[derive(Clone)]
pub struct TensorRSMultilinearPCSConfig<F: FieldExt> {
    pub expansion_factor: usize,
    pub domain_powers: Option<Vec<Vec<F>>>,
    pub fft_domain: Option<Vec<F>>,
    pub ecfft_config: Option<ECFFTConfig<F>>,
    pub l: usize,
    pub num_entries: usize,
    pub num_rows: usize,
}

impl<F: FieldExt> TensorRSMultilinearPCSConfig<F> {
    pub fn num_cols(&self) -> usize {
        self.num_entries / self.num_rows()
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
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
    u_hat_comm: [u8; 32],
    pub test_u_prime: Vec<F>,
    pub test_r_prime: Vec<F>,
    pub eval_r_prime: Vec<F>,
    pub eval_u_prime: Vec<F>,
}

impl<F: FieldExt> TensorMultilinearPCS<F> {
    pub fn new(config: TensorRSMultilinearPCSConfig<F>) -> Self {
        Self { config }
    }

    pub fn commit(&self, poly: &SparseMLPoly<F>) -> CommittedTensorCode<F> {
        // Merkle commit to the evaluations of the polynomial
        let tensor_code = self.encode_zk(poly);
        let tree = tensor_code.commit(self.config.num_cols(), self.config.num_rows());
        tree
    }

    pub fn open(
        &self,
        u_hat_comm: &CommittedTensorCode<F>,
        poly: &SparseMLPoly<F>,
        point: &[F],
        transcript: &mut Transcript<F>,
    ) -> TensorMLOpening<F> {
        let num_cols = self.config.num_cols();
        let num_rows = self.config.num_rows();
        debug_assert_eq!(poly.num_vars, point.len());

        transcript.append_bytes(&u_hat_comm.committed_tree.root());

        // ########################################
        // Testing phase
        // Prove the consistency between the random linear combination of the evaluation tensor (u_prime)
        // and the tensor codeword (u_hat)
        // ########################################

        // Derive the challenge vector;
        let r_u = transcript.challenge_vec(num_rows);

        let u = (0..num_rows)
            .map(|i| {
                poly.evals[(i * num_cols)..((i + 1) * num_cols)]
                    .iter()
                    .map(|entry| entry.1)
                    .collect::<Vec<F>>()
            })
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
        let mut point_rev = point.to_vec();
        point_rev.reverse();

        let log2_num_rows = (num_rows as f64).log2() as usize;
        let q1 = EqPoly::new(point_rev[0..log2_num_rows].to_vec()).evals();

        let eval_r_prime = rlc_rows(blinder, &q1);

        let eval_u_prime = rlc_rows(u.clone(), &q1);

        let eval_queries = self.test_phase(&indices, &u_hat_comm);

        TensorMLOpening {
            x: point.to_vec(),
            y: poly.eval(&point_rev),
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
        }
    }
}

impl<F: FieldExt> TensorMultilinearPCS<F> {
    pub fn verify(
        &self,
        opening: &TensorMLOpening<F>,
        commitment: &[u8; 32],
        transcript: &mut Transcript<F>,
    ) {
        let num_rows = self.config.num_rows();
        let num_cols = self.config.num_cols();

        let u_hat_comm = opening.u_hat_comm;
        transcript.append_bytes(&u_hat_comm);

        assert_eq!(&u_hat_comm, commitment);

        // Verify the base opening

        let base_opening = &opening.base_opening;
        base_opening.verify(u_hat_comm);

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

        let mut x_rev = opening.x.clone();
        x_rev.reverse();

        let log2_num_rows = (num_rows as f64).log2() as usize;
        let q1 = EqPoly::new(x_rev[0..log2_num_rows].to_vec()).evals();
        let q2 = EqPoly::new(x_rev[log2_num_rows..].to_vec()).evals();

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
        let codeword = if self.config.fft_domain.is_some() {
            let fft_domain = self.config.fft_domain.as_ref().unwrap();
            let mut padded_coeffs = message.clone().to_vec();
            padded_coeffs.resize(fft_domain.len(), F::ZERO);
            let codeword = fft(&padded_coeffs, &fft_domain);

            codeword
        } else if self.config.ecfft_config.is_some() {
            let ecfft_config = self.config.ecfft_config.as_ref().unwrap();
            assert_eq!(
                message.len() * self.config.expansion_factor,
                ecfft_config.domain[0].len()
            );
            let extended_evals = extend(
                message,
                &ecfft_config.domain,
                &ecfft_config.matrices,
                &ecfft_config.inverse_matrices,
                0,
            );

            let codeword = [message.to_vec(), extended_evals].concat();
            codeword
        } else {
            let domain_powers = self.config.domain_powers.as_ref().unwrap();
            assert_eq!(message.len(), domain_powers[0].len());
            assert_eq!(
                message.len() * self.config.expansion_factor,
                domain_powers.len()
            );

            let codeword = domain_powers
                .iter()
                .map(|powers| {
                    message
                        .iter()
                        .zip(powers.iter())
                        .fold(F::ZERO, |acc, (m, p)| acc + *m * *p)
                })
                .collect::<Vec<F>>();

            codeword
        };

        codeword
    }

    fn test_phase(&self, indices: &[usize], u_hat_comm: &CommittedTensorCode<F>) -> Vec<Vec<F>> {
        let num_cols = self.config.num_cols() * 2;

        // Query the columns of u_hat
        let num_indices = self.config.l;

        let u_hat_openings = indices
            .iter()
            .map(|index| u_hat_comm.query_column(*index, num_cols))
            .collect::<Vec<Vec<F>>>();

        debug_assert_eq!(u_hat_openings.len(), num_indices);

        u_hat_openings
    }

    fn encode_zk(&self, poly: &SparseMLPoly<F>) -> TensorCode<F> {
        let num_rows = self.config.num_rows();
        let num_cols = self.config.num_cols();

        let codewords = (0..num_rows)
            .map(|i| {
                poly.evals[i * num_cols..(i + 1) * num_cols]
                    .iter()
                    .map(|entry| entry.1)
                    .collect::<Vec<F>>()
            })
            .map(|row| self.split_encode(&row))
            .collect::<Vec<Vec<F>>>();

        TensorCode(codewords)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rs_config::{ecfft, naive, smooth};

    const TEST_NUM_VARS: usize = 10;
    const TEST_L: usize = 10;

    fn test_poly<F: FieldExt>() -> SparseMLPoly<F> {
        let num_entries: usize = 2usize.pow(TEST_NUM_VARS as u32);

        let evals = (0..num_entries)
            .map(|i| (i, F::from(i as u64)))
            .collect::<Vec<(usize, F)>>();

        let ml_poly = SparseMLPoly::new(evals, TEST_NUM_VARS);
        ml_poly
    }

    fn prove_and_verify<F: FieldExt>(ml_poly: SparseMLPoly<F>, pcs: TensorMultilinearPCS<F>) {
        let comm = pcs.commit(&ml_poly);

        let open_at = (0..ml_poly.num_vars)
            .map(|i| F::from(i as u64))
            .collect::<Vec<F>>();

        let mut prover_transcript = Transcript::<F>::new(b"test");
        let opening = pcs.open(&comm, &ml_poly, &open_at, &mut prover_transcript);

        let mut verifier_transcript = Transcript::<F>::new(b"test");
        pcs.verify(
            &opening,
            &comm.committed_tree.root(),
            &mut verifier_transcript,
        );
    }

    fn config_base<F: FieldExt>(ml_poly: &SparseMLPoly<F>) -> TensorRSMultilinearPCSConfig<F> {
        let num_vars = ml_poly.num_vars;
        let num_evals = 2usize.pow(num_vars as u32);
        let num_rows = 2usize.pow((num_vars / 2) as u32);

        let expansion_factor = 2;

        TensorRSMultilinearPCSConfig::<F> {
            expansion_factor,
            domain_powers: None,
            fft_domain: None,
            ecfft_config: None,
            l: TEST_L,
            num_entries: num_evals,
            num_rows,
        }
    }

    #[test]
    fn test_tensor_pcs_fft() {
        type F = halo2curves::pasta::Fp;
        // FFT config
        let ml_poly = test_poly();
        let mut config = config_base(&ml_poly);
        config.fft_domain = Some(smooth::gen_config(config.num_cols()));

        // Test FFT PCS
        let tensor_pcs_fft = TensorMultilinearPCS::<F>::new(config);
        prove_and_verify(ml_poly, tensor_pcs_fft);
    }

    #[test]
    fn test_tensor_pcs_ecfft() {
        type F = halo2curves::secp256k1::Fp;
        let ml_poly = test_poly();

        let mut config = config_base(&ml_poly);
        config.ecfft_config = Some(ecfft::gen_config(config.num_cols()));

        // Test FFT PCS
        let tensor_pcs_ecf = TensorMultilinearPCS::<F>::new(config);
        prove_and_verify(ml_poly, tensor_pcs_ecf);
    }

    #[test]
    fn test_tensor_pcs_naive() {
        type F = halo2curves::secp256k1::Fp;
        // FFT config
        let ml_poly = test_poly();

        // Naive config
        let mut config = config_base(&ml_poly);
        config.domain_powers = Some(naive::gen_config(config.num_cols()));

        // Test FFT PCS
        let tensor_pcs_naive = TensorMultilinearPCS::<F>::new(config);
        prove_and_verify(ml_poly, tensor_pcs_naive);
    }
}
