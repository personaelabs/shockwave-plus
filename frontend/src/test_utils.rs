use crate::constraint_system::ConstraintSystem;
use shockwave_plus::FieldGC;

#[allow(unused_must_use)]
pub fn mock_circuit<F: FieldGC>(num_cons: usize) -> impl Fn(&mut ConstraintSystem<F>) {
    let synthesizer = move |cs: &mut ConstraintSystem<F>| {
        let a = cs.alloc_priv_input();
        let b = cs.alloc_priv_input();

        // There is always one constraint all teh additions,
        // so we do num_cons - 1 multiplications to
        // obtain a circuit with num_cons constraints.
        for _ in 0..(num_cons - 1) {
            a * b;
        }

        cs.expose_public(a * b);
    };

    synthesizer
}

#[allow(unused_must_use)]
#[allow(dead_code)]
pub fn synthetic_circuit<F: FieldGC>() -> (impl Fn(&mut ConstraintSystem<F>), Vec<F>, Vec<F>, Vec<F>)
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
