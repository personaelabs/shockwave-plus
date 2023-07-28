mod fft;
mod polynomial;
pub mod rs_config;
mod tensor_code;
mod tensor_pcs;
mod transcript;
mod tree;
mod utils;
use halo2curves::ff::FromUniformBytes;

pub trait FieldExt: FromUniformBytes<64, Repr = [u8; 32]> {}

impl FieldExt for halo2curves::secp256k1::Fp {}
impl FieldExt for halo2curves::pasta::Fp {}

pub use polynomial::eq_poly::EqPoly;
pub use polynomial::sparse_ml_poly::SparseMLPoly;
pub use tensor_pcs::{TensorMLOpening, TensorMultilinearPCS, TensorRSMultilinearPCSConfig};
pub use transcript::{AppendToTranscript, Transcript};
