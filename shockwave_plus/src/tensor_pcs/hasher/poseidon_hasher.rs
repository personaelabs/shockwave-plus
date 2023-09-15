use super::Hasher;
use crate::{FieldGC, IOPattern, PoseidonSponge};

const PH_WIDTH: usize = 9; // Poseidon hasher width

#[derive(Clone)]
pub struct PoseidonHasher<F: FieldGC> {
    sponge: PoseidonSponge<F, PH_WIDTH>,
}

impl<F: FieldGC> PoseidonHasher<F> {
    pub fn new() -> Self {
        Self {
            sponge: PoseidonSponge::new(b"poseidon_hasher", IOPattern::new(vec![])),
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
