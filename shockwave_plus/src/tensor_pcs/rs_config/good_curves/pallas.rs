use ecfft::GoodCurve;

type Fp = ark_pallas::Fq;

use crate::AppendToTranscript;

use super::FieldGC;

impl AppendToTranscript<Fp> for Fp {
    fn append_to_transcript(&self, transcript: &mut impl crate::TranscriptLike<Fp>) {
        transcript.append_fe(*self);
    }
}

impl FieldGC for Fp {
    fn good_curve(_k: usize) -> (GoodCurve<Fp>, (Fp, Fp)) {
        // Pallas is an FFT-friendly curve, so we don't need an ECFFT curve.
        unimplemented!()
    }
}
