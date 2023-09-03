use ark_ff::{BigInteger, Field, PrimeField};
use merlin::Transcript as MerlinTranscript;
use num_bigint::BigUint;
use std::marker::PhantomData;

#[derive(Clone)]
pub struct Transcript<F: PrimeField> {
    transcript_inner: MerlinTranscript,
    _marker: PhantomData<F>,
}

impl<F: PrimeField> Transcript<F> {
    pub fn new(label: &'static [u8]) -> Self {
        Self {
            transcript_inner: MerlinTranscript::new(label),
            _marker: PhantomData,
        }
    }

    pub fn append_fe(&mut self, fe: &F) {
        self.transcript_inner
            .append_message(b"", &fe.into_bigint().to_bytes_be());
    }

    pub fn append_bytes(&mut self, bytes: &[u8]) {
        self.transcript_inner.append_message(b"", bytes);
    }

    pub fn challenge_vec(&mut self, n: usize) -> Vec<F> {
        (0..n)
            .map(|_| {
                let mut bytes = [0u8; 64];
                self.transcript_inner.challenge_bytes(b"", &mut bytes);
                F::from_random_bytes(&bytes).unwrap()
            })
            .collect()
    }

    pub fn challenge_fe(&mut self) -> F {
        let mut bytes = [0u8; 64];
        self.transcript_inner.challenge_bytes(b"", &mut bytes);
        F::from_random_bytes(&bytes).unwrap()
    }

    pub fn challenge_bytes(&mut self, bytes: &mut [u8]) {
        self.transcript_inner.challenge_bytes(b"", bytes);
    }
}

pub trait AppendToTranscript<F: PrimeField> {
    fn append_to_transcript(&self, transcript: &mut Transcript<F>);
}

impl<F: PrimeField> AppendToTranscript<F> for [u8; 32] {
    fn append_to_transcript(&self, transcript: &mut Transcript<F>) {
        transcript.append_bytes(self);
    }
}
