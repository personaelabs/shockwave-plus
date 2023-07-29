use crate::transcript::Transcript;
use crate::FieldExt;
use tiny_keccak::{Hasher, Keccak};

pub fn rlc_rows<F: FieldExt>(x: Vec<Vec<F>>, r: &[F]) -> Vec<F> {
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

pub fn dot_prod<F: FieldExt>(x: &[F], y: &[F]) -> F {
    assert_eq!(x.len(), y.len());
    let mut result = F::ZERO;
    for i in 0..x.len() {
        result += x[i] * y[i];
    }
    result
}

pub fn hash_two(values: &[[u8; 32]; 2]) -> [u8; 32] {
    let mut hasher = Keccak::v256();
    hasher.update(&values[0]);
    hasher.update(&values[1]);
    let mut hash = [0u8; 32];
    hasher.finalize(&mut hash);
    hash
}

pub fn hash_all(values: &[[u8; 32]]) -> [u8; 32] {
    let mut hasher = Keccak::v256();
    for value in values {
        hasher.update(value);
    }
    let mut hash = [0u8; 32];
    hasher.finalize(&mut hash);
    hash
}

fn sample_index(random_bytes: [u8; 64], size: usize) -> usize {
    let mut acc: u64 = 0;
    for b in random_bytes {
        acc = acc << 8 ^ (b as u64);
    }

    (acc % (size as u64)) as usize
}

pub fn sample_indices<F: FieldExt>(
    num_indices: usize,
    max_index: usize,
    transcript: &mut Transcript<F>,
) -> Vec<usize> {
    assert!(
        num_indices <= max_index,
        "max_index {:?} num_indices {:?}",
        max_index,
        num_indices
    );

    let mut indices = Vec::with_capacity(num_indices);
    let mut counter: u32 = 0;

    let n = max_index / 2;
    while indices.len() < num_indices {
        let mut random_bytes = [0u8; 64];

        transcript.append_bytes(&counter.to_le_bytes());
        transcript.challenge_bytes(&mut random_bytes);

        let index = sample_index(random_bytes, max_index);
        let pair_index = if index > n { index - n } else { index + n };
        if !indices.contains(&index) && !indices.contains(&pair_index) {
            indices.push(index);
        }
        counter += 1;
    }

    indices
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = halo2curves::secp256k1::Fp;

    #[test]
    fn test_sample_indices() {
        let mut transcript = Transcript::<F>::new(b"test_sample_index");
        let num_indices = 10;
        let max_index = 100;
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
