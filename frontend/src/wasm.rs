pub mod wasm_deps {
    // Re-export the dependencies that are used in the wasm module
    pub use crate::constraint_system::{CircuitMeta, ConstraintSystem};
    pub use console_error_panic_hook;
    pub use shockwave_plus::ark_ff::PrimeField;
    pub use shockwave_plus::ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
    pub use shockwave_plus::det_num_cols;
    pub use shockwave_plus::{
        rs_config::{ecfft::gen_config_form_curve, ecfft::ECFFTConfig},
        TensorRSMultilinearPCSConfig,
    };
    pub use shockwave_plus::{
        FieldGC, IOPattern, PoseidonCurve, PoseidonHasher, PoseidonTranscript,
    };

    pub use shockwave_plus::{Proof, ShockwavePlus, R1CS};
    pub use std::sync::Mutex;
    pub use wasm_bindgen;
    pub use wasm_bindgen::prelude::*;

    #[allow(dead_code)]
    pub fn to_felts<F: FieldGC>(bytes: &[u8]) -> Vec<F> {
        bytes
            .chunks_exact(32)
            .map(|x| F::from_be_bytes_mod_order(x))
            .collect::<Vec<F>>()
    }
}

#[allow(unused_imports)]
use wasm_deps::*;

#[macro_export]
macro_rules! circuit {
    ($synthesizer:expr, $field:ty) => {
        const SAMPLE_COLUMNS: usize = 309;
        const EXPANSION_FACTOR: usize = 2;

        static PCS_CONFIG: Mutex<TensorRSMultilinearPCSConfig<$field>> =
            Mutex::new(TensorRSMultilinearPCSConfig {
                expansion_factor: EXPANSION_FACTOR,
                l: SAMPLE_COLUMNS,
                ecfft_config: ECFFTConfig::<$field> {
                    domain: vec![],
                    matrices: vec![],
                    inverse_matrices: vec![],
                },
            });

        static CIRCUIT: Mutex<R1CS<$field>> = Mutex::new(R1CS::empty());

        static CONSTRAINT_SYSTEM: Mutex<ConstraintSystem<$field>> =
            Mutex::new(ConstraintSystem::new());

        pub fn prepare() {
            // ################################
            // Load the circuit
            // ################################

            let mut circuit = CIRCUIT.lock().unwrap();

            let mut cs = CONSTRAINT_SYSTEM.lock().unwrap();
            cs.set_constraints(&$synthesizer);
            *circuit = cs.to_r1cs();

            // ################################
            // Generate the PCS configuration
            // ################################

            let num_cols = det_num_cols(circuit.z_len(), SAMPLE_COLUMNS);
            let k = ((num_cols * EXPANSION_FACTOR).next_power_of_two() as f64).log2() as usize;

            let mut pcs_config = PCS_CONFIG.lock().unwrap();
            let good_curve = <$field>::good_curve(k);

            let ecfft_config = gen_config_form_curve(good_curve.0, good_curve.1);
            pcs_config.ecfft_config = ecfft_config;
        }

        pub fn prove(
            pub_input: &[$field],
            priv_input: &[$field],
        ) -> Proof<$field, PoseidonHasher<$field>> {
            let pcs_config = PCS_CONFIG.lock().unwrap().clone();
            let circuit = CIRCUIT.lock().unwrap().clone();

            let mut cs = CONSTRAINT_SYSTEM.lock().unwrap();
            let witness = cs.gen_witness($synthesizer, pub_input, priv_input);

            // Generate the proof
            let poseidon_hasher = PoseidonHasher::new(PoseidonCurve::SECP256K1);
            let shockwave_plus = ShockwavePlus::new(circuit, pcs_config, poseidon_hasher);

            let mut transcript = PoseidonTranscript::new(
                b"ShockwavePlus",
                PoseidonCurve::SECP256K1,
                IOPattern::new(vec![]),
            );

            let blind = true;
            let proof = shockwave_plus.prove(&witness, pub_input, &mut transcript, blind);
            proof.0
        }

        pub fn verify(proof: Proof<$field, PoseidonHasher<$field>>) -> bool {
            let pcs_config = PCS_CONFIG.lock().unwrap().clone();
            let circuit = CIRCUIT.lock().unwrap().clone();

            let poseidon_hasher = PoseidonHasher::new(PoseidonCurve::SECP256K1);
            let shockwave_plus = ShockwavePlus::new(circuit, pcs_config, poseidon_hasher);

            let mut transcript = PoseidonTranscript::new(
                b"ShockwavePlus",
                PoseidonCurve::SECP256K1,
                IOPattern::new(vec![]),
            );

            shockwave_plus.verify(&proof, &mut transcript);

            true
        }

        // ################################
        // Expose the following functions to the wasm runtime
        // ################################

        #[wasm_bindgen]
        pub fn init_panic_hook() {
            console_error_panic_hook::set_once();
        }

        #[wasm_bindgen]
        pub fn client_prepare() {
            prepare();
        }

        #[wasm_bindgen]
        pub fn client_prove(pub_input: &[u8], priv_input: &[u8]) -> Vec<u8> {
            let pub_input_felts = to_felts(pub_input);
            let priv_input_felts = to_felts(priv_input);

            let proof = prove(&pub_input_felts, &priv_input_felts);
            let mut compressed_bytes = Vec::new();
            proof.serialize_compressed(&mut compressed_bytes).unwrap();
            compressed_bytes
        }

        #[wasm_bindgen]
        pub fn client_verify(proof_ser: &[u8]) -> bool {
            let proof =
                Proof::<$field, PoseidonHasher<$field>>::deserialize_uncompressed_unchecked(
                    proof_ser,
                )
                .unwrap();
            verify(proof)
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::mock_circuit;
    use shockwave_plus::ark_ff::{BigInteger, PrimeField};

    type F = shockwave_plus::ark_secp256k1::Fq;

    #[test]
    fn test_to_felts() {
        let n = 3;
        let felts = (0..n).map(|i| F::from(i)).collect::<Vec<F>>();
        let felt_bytes = felts
            .iter()
            .map(|x| x.into_bigint().to_bytes_be())
            .flatten()
            .collect::<Vec<u8>>();

        let felts_recovered = to_felts::<F>(&felt_bytes);
        assert_eq!(felts, felts_recovered);
    }

    #[test]
    fn test_client_prove() {
        const NUM_CONS: usize = 2usize.pow(8);
        circuit!(mock_circuit(NUM_CONS), F);

        let priv_input = [F::from(3), F::from(4)];
        let pub_input = [priv_input[0] * priv_input[1]];

        prepare();

        let pub_input_bytes = pub_input
            .iter()
            .map(|x| x.into_bigint().to_bytes_be())
            .flatten()
            .collect::<Vec<u8>>();

        let priv_input_bytes = priv_input
            .iter()
            .map(|x| x.into_bigint().to_bytes_be())
            .flatten()
            .collect::<Vec<u8>>();

        let proof_bytes = client_prove(&pub_input_bytes, &priv_input_bytes);

        let proof = Proof::<F, PoseidonHasher<F>>::deserialize_uncompressed_unchecked(
            proof_bytes.as_slice(),
        )
        .unwrap();

        let mut openings_bytes = Vec::new();
        proof
            .z_eval_proof
            .serialize_compressed(&mut openings_bytes)
            .unwrap();
        println!("z_eval_proof.len() = {}", openings_bytes.len());
    }
}
