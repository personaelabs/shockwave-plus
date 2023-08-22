#![allow(non_snake_case)]
mod polynomial;
pub mod rs_config;
mod tensor_code;
mod tensor_rs_pcs;
mod transcript;
mod tree;
mod utils;
use halo2curves::ff::FromUniformBytes;

pub trait FieldExt: FromUniformBytes<64, Repr = [u8; 32]> {}

impl FieldExt for halo2curves::secp256k1::Fp {}
impl FieldExt for halo2curves::pasta::Fp {}

pub use ecfft;
pub use halo2curves;
pub use polynomial::eq_poly::EqPoly;
pub use polynomial::ml_poly::MlPoly;
pub use tensor_rs_pcs::{TensorMLOpening, TensorMultilinearPCS, TensorRSMultilinearPCSConfig};
pub use transcript::{AppendToTranscript, Transcript};
pub use utils::{det_num_cols, det_num_rows, dot_prod};
