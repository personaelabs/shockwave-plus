use crate::{AppendToTranscript, FieldGC, IOPattern, PoseidonCurve, PoseidonSponge};
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

#[derive(Clone)]
pub struct PoseidonHasher<F: FieldGC> {
    sponge: PoseidonSponge<F, 9>,
}

impl<F: FieldGC> PoseidonHasher<F> {
    pub fn new(curve: PoseidonCurve) -> Self {
        Self {
            sponge: PoseidonSponge::new(b"poseidon_hasher", curve, IOPattern::new(vec![])),
        }
    }
}

impl<F: FieldGC> Hasher<F> for PoseidonHasher<F> {
    type T = F;

    fn hash_felts(&self, values: &[F]) -> Self::T {
        let mut sponge = self.sponge.clone();
        sponge.absorb(values);
        sponge.squeeze(1)[0]
    }

    fn hash_all(&self, values: &[Self::T]) -> Self::T {
        self.hash_felts(values)
    }
}
