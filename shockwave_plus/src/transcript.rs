use crate::FieldExt;
use merlin::Transcript as MerlinTranscript;
use std::marker::PhantomData;

#[derive(Clone)]
pub struct Transcript<F: FieldExt> {
    transcript_inner: MerlinTranscript,
    _marker: PhantomData<F>,
}

impl<F: FieldExt> Transcript<F> {
    pub fn new(label: &'static [u8]) -> Self {
        Self {
            transcript_inner: MerlinTranscript::new(label),
            _marker: PhantomData,
        }
    }

    pub fn append_fe(&mut self, fe: &F) {
        self.transcript_inner.append_message(b"", &fe.to_repr());
    }

    pub fn append_bytes(&mut self, bytes: &[u8]) {
        self.transcript_inner.append_message(b"", bytes);
    }

    pub fn challenge_vec(&mut self, n: usize) -> Vec<F> {
        (0..n)
            .map(|_| {
                let mut bytes = [0u8; 64];
                self.transcript_inner.challenge_bytes(b"", &mut bytes);
                F::from_uniform_bytes(&bytes)
            })
            .collect()
    }

    pub fn challenge_fe(&mut self) -> F {
        let mut bytes = [0u8; 64];
        self.transcript_inner.challenge_bytes(b"", &mut bytes);
        F::from_uniform_bytes(&bytes)
    }

    pub fn challenge_bytes(&mut self, bytes: &mut [u8]) {
        self.transcript_inner.challenge_bytes(b"", bytes);
    }
}

pub trait AppendToTranscript<F: FieldExt> {
    fn append_to_transcript(&self, transcript: &mut Transcript<F>);
}

impl<F: FieldExt> AppendToTranscript<F> for [u8; 32] {
    fn append_to_transcript(&self, transcript: &mut Transcript<F>) {
        transcript.append_bytes(self);
    }
}
