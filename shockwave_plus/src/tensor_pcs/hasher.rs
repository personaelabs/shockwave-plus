use crate::{AppendToTranscript, FieldGC, Poseidon, PoseidonCurve};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::fmt::Debug;

pub trait Hasher<F: FieldGC>: Clone + Sync + Send {
    type T: AppendToTranscript<F>
        + Send
        + Sync
        + Clone
        + CanonicalSerialize
        + CanonicalDeserialize
        + PartialEq
        + Debug
        + Copy;
    fn hash_felts(&self, values: &[F]) -> Self::T;
    fn hash_all(&self, values: &[Self::T]) -> Self::T;
}

#[derive(Clone)]
pub struct PoseidonHasher<F: FieldGC> {
    poseidon: Poseidon<F>,
}

impl<F: FieldGC> PoseidonHasher<F> {
    pub fn new(curve: PoseidonCurve) -> Self {
        Self {
            poseidon: Poseidon::new(curve),
        }
    }
}

fn poseidon_hash<F: FieldGC>(input: &[F; 2], poseidon: Poseidon<F>) -> F {
    let mut poseidon = poseidon.clone();
    let h = poseidon.hash(input);
    poseidon.reset();
    h
}

impl<F: FieldGC> Hasher<F> for PoseidonHasher<F> {
    type T = F;

    fn hash_felts(&self, values: &[F]) -> Self::T {
        assert!(values.len().is_power_of_two());

        let m = (values.len() as f64).log2() as usize;
        let mut hashes = values.to_vec();
        let poseidon = self.poseidon.clone();

        for _ in 0..m {
            let n = hashes.len();
            hashes = (0..n)
                .into_par_iter()
                .step_by(2)
                .map(|i| {
                    let h = poseidon_hash(&[hashes[i], hashes[i + 1]], poseidon.clone());
                    h
                })
                .collect::<Vec<F>>();
        }

        // Sanity check
        assert_eq!(hashes.len(), 1);

        hashes[0]
    }

    fn hash_all(&self, values: &[Self::T]) -> Self::T {
        assert!(values.len().is_power_of_two());

        let m = (values.len() as f64).log2() as usize;
        let mut hashes = values.to_vec();
        let poseidon = self.poseidon.clone();

        for _ in 0..m {
            let n = hashes.len();
            hashes = (0..n)
                .into_par_iter()
                .step_by(2)
                .map(|i| {
                    let h = poseidon_hash(&[hashes[i], hashes[i + 1]], poseidon.clone());
                    h
                })
                .collect::<Vec<F>>();
        }

        // Sanity check
        assert_eq!(hashes.len(), 1);

        hashes[0]
    }
}
