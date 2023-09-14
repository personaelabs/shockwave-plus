use blake2b_simd::State;
use std::marker::PhantomData;

use super::Hasher;
use crate::FieldGC;
use ark_ff::BigInteger;

#[derive(Clone)]
pub struct Blake2bHasher<F: FieldGC> {
    _marker: PhantomData<F>,
}

impl<F: FieldGC> Blake2bHasher<F> {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            _marker: PhantomData,
        }
    }
}

impl<F: FieldGC> Hasher<F> for Blake2bHasher<F> {
    type T = [u8; 64];

    fn hash_felts(&self, values: &[F]) -> Self::T {
        let bytes = values
            .iter()
            .map(|v| {
                let mut v_bytes: [u8; 64] = [0; 64];
                v_bytes[..32].copy_from_slice(v.into_bigint().to_bytes_be().as_slice());
                v_bytes
            })
            .collect::<Vec<[u8; 64]>>();
        self.hash_all(&bytes)
    }

    fn hash_all(&self, values: &[Self::T]) -> Self::T {
        let mut state = State::new();

        for val in values {
            state.update(val);
        }
        let result = state.finalize();
        result.as_bytes().try_into().unwrap()
    }
}
