use crate::constraint_system::ConstraintSystem;
use shockwave_plus::ark_ff::PrimeField;

#[allow(unused_must_use)]
pub fn mock_circuit<F: PrimeField>() -> (impl Fn(&mut ConstraintSystem<F>), Vec<F>, Vec<F>, Vec<F>)
{
    let synthesizer = |cs: &mut ConstraintSystem<F>| {
        let w1 = cs.alloc_pub_input();
        let w2 = cs.alloc_pub_input();

        let w6 = cs.alloc_priv_input();

        // Constraint the wires as follows
        // w1 + w2 = w3
        // w1 * w2 = w4
        // w1 * 333 = w5
        // w1 + w6 = w7
        let w3 = w1 + w2;

        w1 * w2;
        cs.mul_const(w1, F::from(333u32));
        cs.add(w1, w6);

        // Expose a wire as a public input
        cs.expose_public(w3);
    };

    let pub_input = vec![F::from(3u32), F::from(4u32), F::from(7u32)];

    let priv_wires_offset = 4;

    // These are the satisfying witnesses values for the above public inputs
    let mut witness = vec![
        F::from(10u32),
        F::from(12u32),
        F::from(999u32),
        F::from(13u32),
    ];
    witness.resize(priv_wires_offset, F::ZERO);

    let priv_input = witness[..1].to_vec();

    (synthesizer, pub_input, priv_input, witness)
}
