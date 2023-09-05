use ark_ff::PrimeField;
use ecfft::GoodCurve;

pub mod secp256k1;

pub trait FieldGC: PrimeField {
    fn good_curve(k: usize) -> (GoodCurve<Self>, (Self, Self));
}
