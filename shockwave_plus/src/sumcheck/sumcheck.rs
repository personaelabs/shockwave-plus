use ark_std::{end_timer, start_timer};

use crate::polynomial::ml_poly::MlPoly;
use crate::sumcheck::unipoly::UniPoly;
use crate::tensor_pcs::{CommittedTensorCode, TensorMultilinearPCS};
use crate::{AppendToTranscript, FieldGC, Hasher, TranscriptLike};

use super::SumCheckProof;

fn init_blinder_poly<F: FieldGC, H: Hasher<F>>(
    num_vars: usize,
    pcs: &TensorMultilinearPCS<F, H>,
    transcript: &mut impl TranscriptLike<F>,
) -> (MlPoly<F>, CommittedTensorCode<F, H>, F) {
    // We implement the zero-knowledge sumcheck protocol
    // described in Section 4.1 https://eprint.iacr.org/2019/317.pdf

    let mut rng = rand::thread_rng();
    // Sample a blinding polynomial g(x_1, ..., x_m)
    let blinder_poly_evals = (0..2usize.pow(num_vars as u32))
        .map(|_| F::rand(&mut rng))
        .collect::<Vec<F>>();
    let blinder_poly = MlPoly::new(blinder_poly_evals.clone());
    let blinder_poly_sum = blinder_poly_evals.iter().fold(F::ZERO, |acc, x| acc + x);

    let commit_b_timer = start_timer!(|| "Commit blinder polynomial");
    let blinder_poly_comm = pcs.commit(&blinder_poly_evals, true);
    end_timer!(commit_b_timer);

    transcript.append_fe(blinder_poly_sum);
    blinder_poly_comm
        .committed_tree
        .root
        .append_to_transcript(transcript);

    (blinder_poly, blinder_poly_comm, blinder_poly_sum)
}

fn challenge_label(label: String) -> String {
    format!("{}-challenge", label)
}

fn rho_label(label: String) -> String {
    format!("{}-rho", label)
}

// This function implements the zero-knowledge sumcheck protocol, and
// is agnostic of the polynomial(s) being summed.
// The function caller must provide the polynomial(s)'s evaluation tables,
// and the function that combines the evaluation tables (i.e. combines the evaluations of polynomials).

// Transcript behavior:
// 1. Append the sum and the commitment to the blinder polynomial.
// 2. Gets a challenge to combine the blinder polynomial with the summed polynomial(s).
// 3. Gets challenges for each round of the sumcheck protocol.

pub fn prove_sum<F: FieldGC, H: Hasher<F>>(
    poly_num_vars: usize,
    poly_degree: usize,
    eval_tables: &mut Vec<Vec<F>>,
    comb_func: impl Fn(&[F]) -> F,
    blind: bool,
    pcs: &TensorMultilinearPCS<F, H>,
    transcript: &mut impl TranscriptLike<F>,
    label: String,
) -> SumCheckProof<F, H> {
    let num_tables = eval_tables.len();
    let mut round_polys = Vec::<UniPoly<F>>::with_capacity(poly_num_vars);

    let (blinder_poly, blinder_poly_comm, blinder_poly_sum) = if blind {
        // We implement the zero-knowledge sumcheck protocol
        // described in Section 4.1 https://eprint.iacr.org/2019/317.pdf
        let (blinder_poly, blinder_poly_comm, blinder_poly_sum) =
            init_blinder_poly(poly_num_vars, pcs, transcript);

        (
            blinder_poly,
            Some(blinder_poly_comm),
            Some(blinder_poly_sum),
        )
    } else {
        (MlPoly::empty(), None, None)
    };

    let mut blinder_table = blinder_poly.evals.clone();

    let rho = if blind {
        Some(transcript.challenge_fe(rho_label(label.clone())))
    } else {
        None
    };

    let challenge = transcript.challenge_vec(poly_num_vars, format!("{}-challenge", label));

    let round_poly_domain = (0..(poly_degree + 1))
        .map(|i| F::from(i as u64))
        .collect::<Vec<F>>();

    let sc_timer = start_timer!(|| "Sumcheck");
    for j in 0..poly_num_vars {
        let high_index = 2usize.pow((poly_num_vars - j - 1) as u32);
        let mut evals = vec![F::ZERO; poly_degree + 1];

        // https://eprint.iacr.org/2019/317.pdf#subsection.3.2
        for b in 0..high_index {
            let r_y_i = challenge[j];
            for (i, eval_at) in round_poly_domain.iter().enumerate() {
                let mut comb_input = Vec::with_capacity(num_tables);
                for table in eval_tables.into_iter() {
                    let table_eval = table[b] + (table[b + high_index] - table[b]) * eval_at;
                    comb_input.push(table_eval);
                }

                evals[i] += comb_func(&comb_input);
            }

            for table in eval_tables.into_iter() {
                table[b] = table[b] + (table[b + high_index] - table[b]) * r_y_i;
            }
        }

        if blind {
            for b in 0..high_index {
                let r_y_i = challenge[j];
                for (i, eval_at) in round_poly_domain.iter().enumerate() {
                    let blinder_eval = blinder_table[b]
                        + (blinder_table[b + high_index] - blinder_table[b]) * eval_at;
                    evals[i] += rho.unwrap() * blinder_eval;
                }

                blinder_table[b] =
                    blinder_table[b] + (blinder_table[b + high_index] - blinder_table[b]) * r_y_i;
            }
        }

        let round_poly = UniPoly::interpolate(&evals);
        round_polys.push(round_poly);
    }

    end_timer!(sc_timer);

    let blinder_poly_eval_proof = if blind {
        Some(pcs.open(
            &blinder_poly_comm.unwrap(),
            &blinder_poly.evals,
            &challenge,
            blinder_poly.eval(&challenge),
            transcript,
            blind,
        ))
    } else {
        None
    };

    SumCheckProof {
        label: label.to_string(),
        blinder_poly_sum,
        round_poly_coeffs: round_polys
            .iter()
            .map(|p| p.coeffs.clone())
            .collect::<Vec<Vec<F>>>(),
        blinder_poly_eval_proof,
    }
}

// Evaluates all the round polynomials at the challenge point,
// and returns the evaluation of the last round polynomial.
pub fn verify_sum<F: FieldGC, H: Hasher<F>>(
    proof: &SumCheckProof<F, H>,
    pcs: &TensorMultilinearPCS<F, H>,
    sum_target: F,
    poly: impl Fn(&[F]) -> F,
    transcript: &mut impl TranscriptLike<F>,
) {
    // Append the sum and the commitment to the blinder polynomial to the transcript.
    if proof.is_blinded() {
        transcript.append_fe(proof.blinder_poly_sum.unwrap());
        proof
            .blinder_poly_eval_proof
            .as_ref()
            .unwrap()
            .u_hat_comm
            .append_to_transcript(transcript);
    }

    // Get the challenge to combine the blinder polynomial with the summed polynomial(s).
    let rho = if proof.is_blinded() {
        Some(transcript.challenge_fe(rho_label(proof.label.clone())))
    } else {
        None
    };

    // Get challenges for each round of the sumcheck protocol.
    let poly_num_vars = proof.round_poly_coeffs.len();
    let challenge = transcript.challenge_vec(poly_num_vars, challenge_label(proof.label.clone()));

    // Verify the validity of the round polynomials.
    let mut target = if proof.is_blinded() {
        sum_target + rho.unwrap() * proof.blinder_poly_sum.unwrap()
    } else {
        sum_target
    };

    for (i, coeffs) in proof.round_poly_coeffs.iter().enumerate() {
        let round_poly = UniPoly::new(coeffs.clone());
        assert_eq!(
            round_poly.eval(F::ZERO) + round_poly.eval(F::ONE),
            target,
            "i = {}",
            i
        );

        target = round_poly.eval(challenge[i]);
    }

    // Verify the opening of the blinder polynomial.
    if proof.is_blinded() {
        pcs.verify(&proof.blinder_poly_eval_proof.as_ref().unwrap(), transcript);
    }

    let poly_eval = if proof.is_blinded() {
        (poly)(&challenge) + rho.unwrap() * proof.blinder_poly_eval_proof.as_ref().unwrap().y
    } else {
        (poly)(&challenge)
    };

    assert_eq!(poly_eval, target);
}

#[cfg(test)]
mod tests {
    use crate::{
        IOPattern, PoseidonCurve, PoseidonHasher, PoseidonTranscript, TensorRSMultilinearPCSConfig,
    };

    use super::*;
    use ark_ff::Field;
    use ark_std::{end_timer, start_timer};
    type Fp = ark_secp256k1::Fq;

    #[test]
    fn test_sumcheck() {
        let poly_num_vars = 8;
        let poly_num_entries = 2usize.pow(poly_num_vars as u32);
        let poly_degree = 2;
        let mut prover_transcript = PoseidonTranscript::new(
            b"test_sumcheck",
            PoseidonCurve::SECP256K1,
            IOPattern::new(vec![]),
        );
        let mut verifier_transcript = prover_transcript.clone();

        let expansion_factor = 2;
        let sample_indices = 309;
        let pcs_config =
            TensorRSMultilinearPCSConfig::new(poly_num_entries, expansion_factor, sample_indices);

        let poseidon_hasher = PoseidonHasher::new(PoseidonCurve::SECP256K1);
        let pcs = TensorMultilinearPCS::new(pcs_config, poseidon_hasher);

        let eval_table_1 = (0..poly_num_entries)
            .map(|i| Fp::from(i as u64))
            .collect::<Vec<Fp>>();

        let eval_table_2 = (0..poly_num_entries)
            .map(|i| Fp::from(i as u64))
            .collect::<Vec<Fp>>();

        let eval_table_3 = eval_table_1
            .iter()
            .zip(eval_table_2.iter())
            .map(|(x, y)| x * y)
            .collect::<Vec<Fp>>();

        let mut eval_tables = vec![
            eval_table_1.clone(),
            eval_table_2.clone(),
            eval_table_3.clone(),
        ];

        let comb_func = |x: &[Fp]| x[0] * x[1] - x[2];

        let blind = false;
        let sumcheck_prove_timer = start_timer!(|| "Sumcheck prove");
        let sumcheck_proof = prove_sum(
            poly_num_vars,
            poly_degree,
            &mut eval_tables,
            comb_func,
            blind,
            &pcs,
            &mut prover_transcript,
            "test_sumcheck".to_string(),
        );
        end_timer!(sumcheck_prove_timer);

        let poly_1 = MlPoly::new(eval_table_1);
        let poly_2 = MlPoly::new(eval_table_2);
        let poly_3 = MlPoly::new(eval_table_3);

        let poly = |x: &[Fp]| poly_1.eval(x) * poly_2.eval(x) - poly_3.eval(x);

        let sum_target = Fp::ZERO;
        verify_sum(
            &sumcheck_proof,
            &pcs,
            sum_target,
            poly,
            &mut verifier_transcript,
        )
    }
}
