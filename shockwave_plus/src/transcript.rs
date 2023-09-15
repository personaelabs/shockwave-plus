use core::panic;
use std::collections::BTreeMap;

use crate::{FieldGC, IOPattern, PoseidonSponge};

pub trait TranscriptLike<F: FieldGC> {
    fn append_fe(&mut self, fe: F);
    fn append_bytes(&mut self, bytes: &[u8]);
    fn challenge_vec(&mut self, n: usize, label: String) -> Vec<F>;
    fn challenge_fe(&mut self, label: String) -> F;
    fn challenge_bytes(&mut self, bytes: &mut [u8], label: String);
    fn get(&self, label: &str) -> F;
}

#[derive(Clone)]
pub struct PoseidonTranscript<F: FieldGC> {
    sponge: PoseidonSponge<F, 3>,
    challenges: BTreeMap<String, F>, // We store challenges for later reference
}

impl<F: FieldGC> PoseidonTranscript<F> {
    pub fn new(label: &'static [u8], io_pattern: IOPattern) -> Self {
        Self {
            sponge: PoseidonSponge::new(label, io_pattern),
            challenges: BTreeMap::new(),
        }
    }
}

impl<F: FieldGC> TranscriptLike<F> for PoseidonTranscript<F> {
    fn append_fe(&mut self, fe: F) {
        self.sponge.absorb(&[fe]);
    }

    fn append_bytes(&mut self, _bytes: &[u8]) {
        let bytes_low = _bytes[0..16].try_into().unwrap();
        let bytes_high = _bytes[16..32].try_into().unwrap();

        let fe_low = F::from_random_bytes(bytes_low).unwrap();
        let fe_high = F::from_random_bytes(bytes_high).unwrap();
        self.sponge.absorb(&[fe_low, fe_high]);
    }

    fn challenge_fe(&mut self, label: String) -> F {
        let c = self.sponge.squeeze(1)[0];
        if label != "".to_string() {
            if self.challenges.contains_key(&label) {
                panic!("Challenge label {} already exists", label);
            }
            self.challenges.insert(label, c);
        }

        c
    }

    fn challenge_bytes(&mut self, _bytes: &mut [u8], _label: String) {
        unimplemented!()
    }

    fn challenge_vec(&mut self, n: usize, label: String) -> Vec<F> {
        let c = self.sponge.squeeze(n);

        for i in 0..n {
            let label_i = format!("{}-{}", label, i);
            if self.challenges.contains_key(label_i.as_str()) {
                panic!("Challenge label {} already exists", label_i);
            }
            self.challenges.insert(label_i, c[i]);
        }
        c
    }

    fn get(&self, label: &str) -> F {
        *self
            .challenges
            .get(label)
            .unwrap_or_else(|| panic!("Challenge label {} does not exist", label))
    }
}

pub trait AppendToTranscript<F: FieldGC> {
    fn append_to_transcript(&self, transcript: &mut impl TranscriptLike<F>);
}

impl<F: FieldGC> AppendToTranscript<F> for [u8; 32] {
    fn append_to_transcript(&self, transcript: &mut impl TranscriptLike<F>) {
        transcript.append_bytes(self);
    }
}

impl<F: FieldGC> AppendToTranscript<F> for [u8; 64] {
    fn append_to_transcript(&self, transcript: &mut impl TranscriptLike<F>) {
        transcript.append_bytes(self[..32].try_into().unwrap());
        transcript.append_bytes(self[32..].try_into().unwrap());
    }
}
