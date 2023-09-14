mod blake2b_hasher;
mod keccak_hasher;
mod poseidon_hasher;
pub use blake2b_hasher::Blake2bHasher;
pub use keccak_hasher::KeccakHasher;
pub use poseidon_hasher::PoseidonHasher;

use crate::{AppendToTranscript, FieldGC};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
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
