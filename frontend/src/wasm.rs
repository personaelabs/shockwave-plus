// Re-export the dependencies that are used in the wasm module
pub use crate::constraint_system::ConstraintSystem;
pub use bincode;
pub use console_error_panic_hook;
pub use shockwave_plus::tensor_pcs::{
    rs_config::{
        ecfft::gen_config_form_curve, ecfft::ECFFTConfig,
        good_curves::secp256k1::secp256k1_good_curve,
    },
    TensorRSMultilinearPCSConfig,
};
pub use shockwave_plus::{FieldExt, Transcript};
pub use shockwave_plus::{PartialSpartanProof, ShockwavePlus, R1CS};
pub use std::sync::Mutex;
pub use wasm_bindgen;
pub use wasm_bindgen::prelude::*;

#[allow(dead_code)]
pub fn to_felts<F: FieldExt>(bytes: &[u8]) -> Vec<F> {
    bytes
        .chunks_exact(32)
        .map(|x| F::from_repr(x.try_into().unwrap()).unwrap())
        .collect::<Vec<F>>()
}

#[macro_export]
macro_rules! test_circuit {
    ($synthesizer:expr, $field:ty) => {
        type F = $field;

        pub fn mock_run(pub_input: &[F], priv_input: &[F]) {
            let mut cs = ConstraintSystem::new();
            let witness = cs.gen_witness($synthesizer, pub_input, priv_input);

            let z = R1CS::construct_z(&witness, pub_input);
            assert!(cs.is_sat(&z, $synthesizer));
        }
    };
}

#[macro_export]
macro_rules! circuit {
    ($synthesizer:expr, $field:ty) => {
        type F = $field;

        static PCS_CONFIG: Mutex<TensorRSMultilinearPCSConfig<F>> =
            Mutex::new(TensorRSMultilinearPCSConfig {
                expansion_factor: 2,
                l: 1,
                ecfft_config: ECFFTConfig::<F> {
                    domain: vec![],
                    matrices: vec![],
                    inverse_matrices: vec![],
                },
            });

        static CIRCUIT: Mutex<R1CS<F>> = Mutex::new(R1CS::empty());

        static CONSTRAINT_SYSTEM: Mutex<ConstraintSystem<F>> = Mutex::new(ConstraintSystem::new());

        pub fn prepare() {
            // ################################
            // Generate the PCS configuration
            // ################################

            let mut pcs_config = PCS_CONFIG.lock().unwrap();
            let good_curve = secp256k1_good_curve(10);

            let ecfft_config = gen_config_form_curve(good_curve.0, good_curve.1);
            pcs_config.ecfft_config = ecfft_config;

            // ################################
            // Load the circuit
            // ################################

            let mut circuit = CIRCUIT.lock().unwrap();

            let mut cs = CONSTRAINT_SYSTEM.lock().unwrap();
            *circuit = cs.gen_constraints($synthesizer);
        }

        pub fn prove(pub_input: &[F], priv_input: &[F]) -> PartialSpartanProof<F> {
            let pcs_config = PCS_CONFIG.lock().unwrap().clone();
            let circuit = CIRCUIT.lock().unwrap().clone();

            let mut cs = CONSTRAINT_SYSTEM.lock().unwrap();
            let witness = cs.gen_witness($synthesizer, pub_input, priv_input);

            // Generate the proof
            let shockwave_plus = ShockwavePlus::new(circuit, pcs_config);

            let mut transcript = Transcript::new(b"ShockwavePlus");
            let proof = shockwave_plus.prove(&witness, pub_input, &mut transcript);
            proof.0
        }

        pub fn verify(proof: PartialSpartanProof<F>) -> bool {
            let pcs_config = PCS_CONFIG.lock().unwrap().clone();
            let circuit = CIRCUIT.lock().unwrap().clone();

            let shockwave_plus = ShockwavePlus::new(circuit, pcs_config);
            let mut verifier_transcript = Transcript::new(b"ShockwavePlus");
            shockwave_plus.verify_partial(&proof, &mut verifier_transcript);

            true
        }

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
            bincode::serialize(&proof).unwrap()
        }

        #[wasm_bindgen]
        pub fn client_verify(proof: &[u8]) -> bool {
            let proof: PartialSpartanProof<F> = bincode::deserialize(proof).unwrap();
            verify(proof)
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::mock_circuit;
    use shockwave_plus::halo2curves::ff::PrimeField;
    use shockwave_plus::halo2curves::secp256k1::Fp;

    type F = Fp;

    #[test]
    fn test_to_felts() {
        let n = 3;
        let felts = (0..n).map(|i| F::from(i)).collect::<Vec<F>>();
        let felt_bytes = felts
            .iter()
            .map(|x| x.to_repr())
            .flatten()
            .collect::<Vec<u8>>();

        let felts_recovered = to_felts::<F>(&felt_bytes);
        assert_eq!(felts, felts_recovered);
    }

    #[test]
    fn test_client_prove() {
        let (_, pub_input, priv_input, witness) = mock_circuit::<F>();
        circuit!(
            |cs: &mut ConstraintSystem<F>| {
                let a_w1 = cs.alloc_pub_input();
                let a_w2 = cs.alloc_pub_input();

                let a_w6 = cs.alloc_priv_input();

                // Constraint the wires as follows
                // w1 + w2 = w3
                // w1 * w2 = w4
                // w1 * 333 = w5
                // w1 + w6 = w7
                let a_w3 = cs.add(a_w1, a_w2);
                let a_w4 = cs.mul(a_w1, a_w2);
                let a_w5 = cs.mul_const(a_w1, F::from(333));
                let a_w7 = cs.add(a_w1, a_w6);

                // Expose the wires as public inputs
                cs.expose_public(a_w3);
            },
            Fp
        );

        prepare();

        let pub_input_bytes = pub_input
            .iter()
            .map(|x| x.to_repr().to_vec())
            .flatten()
            .collect::<Vec<u8>>();

        let priv_input_bytes = priv_input
            .iter()
            .map(|x| x.to_repr().to_vec())
            .flatten()
            .collect::<Vec<u8>>();

        let proof_bytes = client_prove(&pub_input_bytes, &priv_input_bytes);

        let partial_proof: PartialSpartanProof<F> =
            bincode::deserialize(proof_bytes.as_slice()).unwrap();

        let shockwave_plus = ShockwavePlus::new(
            CIRCUIT.lock().unwrap().clone(),
            PCS_CONFIG.lock().unwrap().clone(),
        );

        let mut verifier_transcript = Transcript::new(b"ShockwavePlus");
        shockwave_plus.verify_partial(&partial_proof, &mut verifier_transcript);
    }
}
