use std::marker::PhantomData;

use super::hasher::Hasher;
use crate::rs_config::ecfft::{gen_config_form_curve, ECFFTConfig};
use crate::{timer_end, timer_start, FieldGC};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ecfft::extend;
use rand::thread_rng;
#[cfg(feature = "parallel")]
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use super::tensor_code::TensorCode;
use super::utils::{det_num_cols, det_num_rows, dot_prod, rlc_rows, sample_indices};
use crate::polynomial::eq_poly::EqPoly;
use crate::transcript::TranscriptLike;

use super::tensor_code::CommittedTensorCode;

#[derive(Clone)]
pub struct TensorRSMultilinearPCSConfig<F: FieldGC> {
    pub expansion_factor: usize,
    pub ecfft_config: ECFFTConfig<F>,
    pub l: usize,
}

impl<F: FieldGC> TensorRSMultilinearPCSConfig<F> {
    pub fn new(num_entries: usize, expansion_factor: usize, l: usize) -> Self {
        let num_cols = det_num_cols(num_entries, l);
        let k = ((num_cols * expansion_factor).next_power_of_two() as f64).log2() as usize;

        let (good_curve, coset_offset) = F::good_curve(k);
        let ecfft_config = gen_config_form_curve(good_curve, coset_offset);

        Self {
            expansion_factor,
            ecfft_config,
            l,
        }
    }

    pub fn num_cols(&self, num_entries: usize) -> usize {
        det_num_cols(num_entries, self.l)
    }

    pub fn num_rows(&self, num_entries: usize) -> usize {
        det_num_rows(num_entries, self.l)
    }
}

#[derive(Clone)]
pub struct TensorMultilinearPCS<F: FieldGC, H: Hasher<F>> {
    config: TensorRSMultilinearPCSConfig<F>,
    hasher: H,
    _marker: PhantomData<H>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct TensorMLOpening<F: FieldGC, H: Hasher<F>> {
    pub x: Vec<F>,
    pub y: F,
    pub column_roots: Vec<H::T>,
    pub query_indices: Vec<usize>,
    pub eval_query_leaves: Vec<Vec<F>>,
    pub u_hat_comm: H::T,
    pub eval_r_prime: Option<Vec<F>>,
    pub eval_u_prime: Vec<F>,
    pub poly_num_vars: usize,
}

impl<F: FieldGC, H: Hasher<F>> TensorMLOpening<F, H> {
    pub fn is_blinded(&self) -> bool {
        self.eval_r_prime.is_some()
    }
}

impl<F: FieldGC, H: Hasher<F>> TensorMultilinearPCS<F, H> {
    pub fn new(config: TensorRSMultilinearPCSConfig<F>, hasher: H) -> Self {
        Self {
            config,
            hasher,
            _marker: PhantomData,
        }
    }

    pub fn commit(&self, ml_poly_evals: &[F], blind: bool) -> CommittedTensorCode<F, H> {
        let n = ml_poly_evals.len();
        assert!(n.is_power_of_two());

        let encode_timer = timer_start("Encoding");
        let tensor_code = self.encode(ml_poly_evals, blind);
        timer_end(encode_timer);
        let tree = tensor_code.commit(
            self.config.num_cols(n),
            self.config.num_rows(n),
            &self.hasher,
        );
        tree
    }

    pub fn open(
        &self,
        u_hat_comm: &CommittedTensorCode<F, H>,
        // TODO: Remove poly and use u_hat_comm
        ml_poly_evals: &[F],
        point: &[F],
        eval: F,
        transcript: &mut impl TranscriptLike<F>,
        blind: bool,
    ) -> TensorMLOpening<F, H> {
        let n = ml_poly_evals.len();
        assert!(n.is_power_of_two());
        let num_vars = (n as f64).log2() as usize;

        let num_cols = self.config.num_cols(n);
        let num_rows = self.config.num_rows(n);

        debug_assert_eq!(num_vars, point.len());

        let u = (0..num_rows)
            .map(|i| ml_poly_evals[(i * num_cols)..((i + 1) * num_cols)].to_vec())
            .collect::<Vec<Vec<F>>>();

        let num_indices = self.config.l;
        let indices = sample_indices(num_indices, num_cols * 2, transcript);

        // ########################################
        // Evaluation phase
        // Prove the consistency
        // ########################################

        // Get the evaluation point

        let log2_num_rows = (num_rows as f64).log2() as usize;
        let q1 = EqPoly::new(point[0..log2_num_rows].to_vec()).evals();

        let eval_r_prime = if blind {
            // Random linear combination of the blinder
            let blinder = u_hat_comm
                .tensor_codeword
                .0
                .iter()
                .map(|row| row[(row.len() / 2)..].to_vec())
                .collect::<Vec<Vec<F>>>();

            Some(rlc_rows(blinder, &q1))
        } else {
            None
        };

        let eval_u_prime = rlc_rows(u.clone(), &q1);

        let eval_queries = self.test_phase(&indices, &u_hat_comm);

        TensorMLOpening {
            x: point.to_vec(),
            y: eval,
            query_indices: indices,
            eval_query_leaves: eval_queries,
            u_hat_comm: u_hat_comm.committed_tree.root(),
            eval_r_prime,
            eval_u_prime,
            column_roots: u_hat_comm.committed_tree.column_roots.clone(),
            poly_num_vars: num_vars,
        }
    }

    pub fn verify(&self, opening: &TensorMLOpening<F, H>, transcript: &mut impl TranscriptLike<F>) {
        let poly_num_entries = 2usize.pow(opening.poly_num_vars as u32);
        let num_rows = self.config.num_rows(poly_num_entries);
        let num_cols = self.config.num_cols(poly_num_entries);

        // Verify that the column roots hashes to the commit to the tensor code
        assert_eq!(
            opening.u_hat_comm,
            self.hasher.hash_all(&opening.column_roots)
        );

        // ########################################
        // Verify evaluation phase
        // ########################################

        let num_indices = self.config.l;
        let indices = sample_indices(num_indices, num_cols * 2, transcript);

        let log2_num_rows = (num_rows as f64).log2() as usize;
        let q1 = EqPoly::new(opening.x[0..log2_num_rows].to_vec()).evals();
        let q2 = EqPoly::new(opening.x[log2_num_rows..].to_vec()).evals();

        let mut eval_u_prime_rs_codeword = self.rs_encode(&opening.eval_u_prime);

        if opening.is_blinded() {
            // The proof is zk
            eval_u_prime_rs_codeword = eval_u_prime_rs_codeword
                .iter()
                .zip(opening.eval_r_prime.as_ref().unwrap().iter())
                .map(|(c, r)| *c + *r)
                .collect::<Vec<F>>();
        }

        debug_assert_eq!(q1.len(), opening.eval_query_leaves[0].len());

        #[cfg(feature = "parallel")]
        indices
            .par_iter()
            .zip(opening.eval_query_leaves.par_iter())
            .for_each(|(expected_index, leaves)| {
                // Verify that the opened column hashes to the column root
                let column_root = self.hasher.hash_felts(&leaves);
                let expected_column_root = opening.column_roots[*expected_index];
                assert_eq!(column_root, expected_column_root);

                let mut sum = F::ZERO;
                for (leaf, q1_i) in leaves.iter().zip(q1.iter()) {
                    sum += *q1_i * *leaf;
                }
                assert_eq!(sum, eval_u_prime_rs_codeword[*expected_index]);
            });

        #[cfg(not(feature = "parallel"))]
        indices
            .iter()
            .zip(opening.eval_query_leaves.iter())
            .for_each(|(expected_index, leaves)| {
                // Verify that the opened column hashes to the column root
                let column_root = self.hasher.hash_felts(&leaves);
                let expected_column_root = opening.column_roots[*expected_index];
                assert_eq!(column_root, expected_column_root);

                let mut sum = F::ZERO;
                for (leaf, q1_i) in leaves.iter().zip(q1.iter()) {
                    sum += *q1_i * *leaf;
                }
                assert_eq!(sum, eval_u_prime_rs_codeword[*expected_index]);
            });

        let expected_eval = dot_prod(&opening.eval_u_prime, &q2);
        assert_eq!(expected_eval, opening.y);
    }

    fn split_encode(&self, message: &[F]) -> Vec<F> {
        let codeword = self.rs_encode(message);

        let mut rng = thread_rng();
        let blinder = (0..codeword.len())
            .map(|_| F::rand(&mut rng))
            .collect::<Vec<F>>();

        let mut randomized_codeword = codeword
            .iter()
            .zip(blinder.iter())
            .map(|(c, b)| *b + *c)
            .collect::<Vec<F>>();

        randomized_codeword.extend_from_slice(&blinder);
        debug_assert_eq!(randomized_codeword.len(), codeword.len() * 2);
        randomized_codeword
    }

    fn rs_encode(&self, message: &[F]) -> Vec<F> {
        assert!(message.len().is_power_of_two());

        let codeword_len = message.len() * self.config.expansion_factor;

        let codeword_len_log2 = (codeword_len as f64).log2() as usize;

        // Resize the domain to the correct size

        let ecfft_config = &self.config.ecfft_config;
        let config_domain_size = ecfft_config.domain.len();

        debug_assert!(config_domain_size >= codeword_len_log2 - 1);
        let matrices =
            ecfft_config.matrices[(config_domain_size - (codeword_len_log2 - 1))..].to_vec();
        let inverse_matrices = ecfft_config.inverse_matrices
            [(config_domain_size - (codeword_len_log2 - 1))..]
            .to_vec();

        let extended_evals = extend(&message, &matrices, &inverse_matrices, 0);

        let codeword = [message.to_vec(), extended_evals].concat();

        codeword
    }

    fn test_phase(&self, indices: &[usize], u_hat_comm: &CommittedTensorCode<F, H>) -> Vec<Vec<F>> {
        // Query the columns of u_hat
        let num_indices = self.config.l;

        let u_hat_openings = indices
            .iter()
            .map(|index| u_hat_comm.query_column(*index))
            .collect::<Vec<Vec<F>>>();

        debug_assert_eq!(u_hat_openings.len(), num_indices);

        u_hat_openings
    }

    fn encode(&self, ml_poly_evals: &[F], blind: bool) -> TensorCode<F> {
        let n = ml_poly_evals.len();
        assert!(n.is_power_of_two());

        let num_rows = self.config.num_rows(n);
        let num_cols = self.config.num_cols(n);
        debug_assert_eq!(n, num_cols * num_rows);

        let rows = (0..num_rows)
            .map(|i| &ml_poly_evals[i * num_cols..(i + 1) * num_cols])
            .collect::<Vec<&[F]>>();

        let codewords = if blind {
            #[cfg(feature = "parallel")]
            let codewords = rows
                .par_iter()
                .map(|row| self.split_encode(row))
                .collect::<Vec<Vec<F>>>();

            #[cfg(not(feature = "parallel"))]
            let codewords = rows
                .iter()
                .map(|row| self.split_encode(row))
                .collect::<Vec<Vec<F>>>();

            codewords
        } else {
            #[cfg(feature = "parallel")]
            let codewords = rows
                .par_iter()
                .map(|row| self.rs_encode(row))
                .collect::<Vec<Vec<F>>>();

            #[cfg(not(feature = "parallel"))]
            let codewords = rows
                .iter()
                .map(|row| self.rs_encode(row))
                .collect::<Vec<Vec<F>>>();

            codewords
        };

        TensorCode(codewords)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        polynomial::ml_poly::MlPoly, tensor_pcs::hasher::PoseidonHasher,
        transcript::PoseidonTranscript, AppendToTranscript, IOPattern,
    };

    const TEST_NUM_VARS: usize = 8;
    const TEST_L: usize = 10;

    fn test_poly_evals<F: FieldGC>() -> MlPoly<F> {
        let num_entries: usize = 2usize.pow(TEST_NUM_VARS as u32);

        let evals = (0..num_entries)
            .map(|i| F::from((i + 1) as u64))
            .collect::<Vec<F>>();

        MlPoly::new(evals)
    }

    fn prove_and_verify<F: FieldGC, H: Hasher<F>>(
        ml_poly: &MlPoly<F>,
        pcs: TensorMultilinearPCS<F, H>,
    ) {
        let ml_poly_evals = &ml_poly.evals;
        let blind = true;
        let comm = pcs.commit(ml_poly_evals, blind);

        let ml_poly_num_vars = (ml_poly_evals.len() as f64).log2() as usize;
        let open_at = (0..ml_poly_num_vars)
            .map(|i| F::from(i as u64))
            .collect::<Vec<F>>();
        let y = ml_poly.eval(&open_at);

        let mut prover_transcript = PoseidonTranscript::new(b"test", IOPattern::new(vec![]));

        comm.committed_tree
            .root
            .append_to_transcript(&mut prover_transcript);

        let opening = pcs.open(
            &comm,
            ml_poly_evals,
            &open_at,
            y,
            &mut prover_transcript,
            blind,
        );

        let mut verifier_transcript = PoseidonTranscript::new(b"test", IOPattern::new(vec![]));
        comm.committed_tree
            .root
            .append_to_transcript(&mut verifier_transcript);

        pcs.verify(&opening, &mut verifier_transcript);
    }

    #[test]
    fn test_tensor_pcs() {
        type F = ark_secp256k1::Fq;
        let ml_poly = test_poly_evals();

        let expansion_factor = 2;

        let n = ml_poly.evals.len();
        let config = TensorRSMultilinearPCSConfig::<F>::new(n, expansion_factor, TEST_L);

        let poseidon_hasher = PoseidonHasher::new();
        let tensor_pcs_ecf = TensorMultilinearPCS::<F, _>::new(config, poseidon_hasher);
        prove_and_verify(&ml_poly, tensor_pcs_ecf);
    }
}
