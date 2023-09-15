use crate::transcript::TranscriptLike;
use crate::FieldGC;
use ark_ff::BigInteger;

pub fn rlc_rows<F: FieldGC>(x: Vec<Vec<F>>, r: &[F]) -> Vec<F> {
    debug_assert_eq!(x.len(), r.len());
    let num_cols = x[0].len();
    let mut result = vec![F::ZERO; num_cols];
    for (row, r_i) in x.iter().zip(r.iter()) {
        for j in 0..num_cols {
            result[j] += row[j] * r_i
        }
    }
    result
}

pub fn dot_prod<F: FieldGC>(x: &[F], y: &[F]) -> F {
    debug_assert_eq!(x.len(), y.len());
    let mut result = F::ZERO;
    for i in 0..x.len() {
        result += x[i] * y[i];
    }
    result
}

fn sample_index<F: FieldGC>(transcript: &mut impl TranscriptLike<F>, max_index: usize) -> usize {
    let max_index_log2 = (max_index as f64).log2() as usize;
    let challenge = transcript.challenge_fe("".to_string());
    let challenge_bits = challenge.into_bigint().to_bits_be();

    let bits = &challenge_bits[..max_index_log2];

    let mut acc: usize = 0;
    for b in bits {
        acc <<= 1;
        acc += *b as usize;
    }

    acc
}

pub fn sample_indices<F: FieldGC>(
    num_indices: usize,
    max_index: usize,
    transcript: &mut impl TranscriptLike<F>,
) -> Vec<usize> {
    assert!(max_index.is_power_of_two());

    assert!(
        num_indices <= max_index,
        "max_index {:?} num_indices {:?}",
        max_index,
        num_indices
    );

    let mut indices = Vec::with_capacity(num_indices);
    let mut counter = F::ZERO;

    let n = max_index / 2;
    while indices.len() < num_indices {
        transcript.append_fe(counter);

        let index = sample_index(transcript, max_index);
        let pair_index = if index > n { index - n } else { index + n };
        if !indices.contains(&index) && !indices.contains(&pair_index) {
            indices.push(index);
        }
        counter += F::ONE;
    }

    indices
}

pub fn det_num_cols(num_entries: usize, l: usize) -> usize {
    assert!(num_entries.is_power_of_two());
    let num_entries_sqrt = (num_entries as f64).sqrt() as usize;
    // The number of columns must be a power of two
    // to tensor-query the polynomial evaluation
    let num_cols = std::cmp::max(num_entries_sqrt, l).next_power_of_two();
    num_cols
}

pub fn det_num_rows(num_entries: usize, l: usize) -> usize {
    assert!(num_entries.is_power_of_two());
    // The number of rows must be a power of two
    // to tensor-query the polynomial evaluation
    let num_rows = (num_entries / det_num_cols(num_entries, l)).next_power_of_two();
    num_rows
}

#[cfg(test)]
mod tests {
    use crate::{transcript::PoseidonTranscript, IOPattern};

    use super::*;
    type F = ark_secp256k1::Fq;

    #[test]
    fn test_sample_indices() {
        let mut transcript =
            PoseidonTranscript::<F>::new(b"test_sample_index", IOPattern::new(vec![]));
        let num_indices = 10;
        let max_index = 128;
        let indices = sample_indices(num_indices, max_index, &mut transcript);

        assert_eq!(indices.len(), 10);
        let n = max_index / 2;
        for index in &indices {
            if *index > n {
                assert!(!indices.contains(&(index - n)));
            } else {
                assert!(!indices.contains(&(index + n)));
            }
        }
    }
}
