use ark_ff::PrimeField;
use ecfft::GoodCurve;

use crate::AppendToTranscript;

pub mod secp256k1;

pub trait FieldGC: PrimeField + Sync + Send + AppendToTranscript<Self> {
    fn good_curve(k: usize) -> (GoodCurve<Self>, (Self, Self));
}
