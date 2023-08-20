use crate::{constraint_system::ConstraintSystem, FieldExt};

pub fn mock_circuit<F: FieldExt>() -> (impl Fn(&mut ConstraintSystem<F>), Vec<F>, Vec<F>, Vec<F>) {
    let synthesizer = |cs: &mut ConstraintSystem<F>| {
        let a_w1 = cs.alloc_pub_input("w1");
        let a_w2 = cs.alloc_pub_input("w2");

        let a_w6 = cs.alloc_priv_input("w6");

        // Constraint the wires as follows
        // w1 + w2 = w3
        // w1 * w2 = w4
        // w1 * 333 = w5
        // w1 + w6 = w7
        cs.add(a_w1, a_w2, "w3");
        cs.mul(a_w1, a_w2, "w4");
        cs.mul_const(a_w1, F::from(333), "w5");
        cs.add(a_w1, a_w6, "w7");

        // Expose a wire as a public input
        cs.expose_public("w3");
    };

    let pub_input = vec![F::from(3), F::from(4), F::from(7)];

    let priv_wires_offset = 4;

    // These are the satisfying witnesses values for the above public inputs
    let mut witness = vec![F::from(10), F::from(12), F::from(999), F::from(13)];
    witness.resize(priv_wires_offset, F::ZERO);

    let priv_input = witness[..1].to_vec();

    (synthesizer, pub_input, priv_input, witness)
}