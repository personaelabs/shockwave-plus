use crate::{FieldGC, IOPattern, PoseidonCurve, PoseidonSponge};

pub trait TranscriptLike<F: FieldGC> {
    fn append_fe(&mut self, fe: F);
    fn append_bytes(&mut self, bytes: &[u8]);
    fn challenge_vec(&mut self, n: usize) -> Vec<F>;
    fn challenge_fe(&mut self) -> F;
    fn challenge_bytes(&mut self, bytes: &mut [u8]);
}

pub struct PoseidonTranscript<F: FieldGC> {
    sponge: PoseidonSponge<F>,
}

impl<F: FieldGC> PoseidonTranscript<F> {
    pub fn new(label: &'static [u8], curve: PoseidonCurve, io_pattern: IOPattern) -> Self {
        Self {
            sponge: PoseidonSponge::new(label, curve, io_pattern),
        }
    }
}

impl<F: FieldGC> TranscriptLike<F> for PoseidonTranscript<F> {
    fn append_fe(&mut self, fe: F) {
        self.sponge.absorb(&[fe]);
    }

    fn append_bytes(&mut self, _bytes: &[u8]) {
        unimplemented!()
    }

    fn challenge_fe(&mut self) -> F {
        self.sponge.squeeze(1)[0]
    }

    fn challenge_bytes(&mut self, _bytes: &mut [u8]) {
        unimplemented!()
    }

    fn challenge_vec(&mut self, n: usize) -> Vec<F> {
        self.sponge.squeeze(n)
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
