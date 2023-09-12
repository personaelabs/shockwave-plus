use std::marker::PhantomData;

use tiny_keccak::{Hasher as TKeccakHasher, Keccak};

use super::Hasher;
use crate::FieldGC;
use ark_ff::BigInteger;

#[derive(Clone)]
pub struct KeccakHasher<F: FieldGC> {
    _marker: PhantomData<F>,
}

impl<F: FieldGC> KeccakHasher<F> {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            _marker: PhantomData,
        }
    }
}

impl<F: FieldGC> Hasher<F> for KeccakHasher<F> {
    type T = [u8; 32];

    fn hash_felts(&self, values: &[F]) -> Self::T {
        let bytes = values
            .iter()
            .map(|v| v.into_bigint().to_bytes_be().try_into().unwrap())
            .collect::<Vec<[u8; 32]>>();
        self.hash_all(&bytes)
    }

    fn hash_all(&self, values: &[Self::T]) -> Self::T {
        let mut hasher = Keccak::v256();
        for val in values {
            hasher.update(val);
        }
        let mut result = [0u8; 32];
        hasher.finalize(&mut result);
        result
    }
}
