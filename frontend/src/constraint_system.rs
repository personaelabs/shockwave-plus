use std::marker::PhantomData;

use shockwave_plus::{Matrix, SparseMatrixEntry, R1CS};

use crate::FieldExt;

#[derive(Debug, Clone, Copy)]
pub struct Wire<F: FieldExt> {
    id: usize,
    index: usize,
    label: &'static str,
    _marker: PhantomData<F>,
}

impl<F: FieldExt> Wire<F> {
    pub fn new(id: usize, index: usize) -> Self {
        Wire {
            id,
            index,
            label: "",
            _marker: PhantomData,
        }
    }

    pub fn assign_label(&mut self, label: &'static str) {
        self.label = label;
    }

    pub fn label(&self) -> &'static str {
        if self.label == "" {
            "no label"
        } else {
            self.label
        }
    }

    pub fn add(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.add(*self, w)
    }

    pub fn add_const(&self, c: F, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.add_const(*self, c)
    }

    pub fn neg(&self, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.neg(*self)
    }

    pub fn sub(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.sub(*self, w)
    }

    pub fn sub_const(&self, c: F, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.sub_const(*self, c)
    }

    pub fn mul(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.mul(*self, w)
    }

    pub fn mul_const(&self, c: F, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.mul_const(*self, c)
    }

    pub fn square(&self, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.square(*self)
    }

    pub fn div(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.div(*self, w)
    }

    pub fn div_or_zero(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.div_or_zero(*self, w)
    }

    pub fn is_zero(&self, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.is_zero(*self)
    }

    pub fn is_equal(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.is_equal(*self, w)
    }

    pub fn assert_zero(&self, cs: &mut ConstraintSystem<F>) {
        cs.assert_zero(*self)
    }

    pub fn assert_equal(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) {
        cs.assert_equal(*self, w)
    }

    pub fn not(&self, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.not(*self)
    }

    pub fn val(&self, cs: &mut ConstraintSystem<F>) -> Option<F> {
        if cs.mode == Mode::WitnessGen {
            Some(cs.wires[self.index])
        } else {
            None
        }
    }
}

// Stores a linear combination a_0 * w_0 + a_1 * w_1 + ... + a_{n-1} * w_{n-1}
// in a sparse vector of tuples.
#[derive(Debug, Clone)]
pub struct LinearCombination<F>(Vec<F>);

impl<F: FieldExt> LinearCombination<F> {
    pub fn new(max_terms: usize) -> Self {
        let terms = vec![F::ZERO; max_terms];
        Self(terms)
    }

    pub fn increment_coeff(&mut self, wire: Wire<F>) {
        self.0[wire.index] += F::ONE;
    }

    pub fn increment_coeff_by(&mut self, wire: Wire<F>, inc: F) {
        self.0[wire.index] += inc;
    }

    pub fn set_coeff(&mut self, wire: Wire<F>, coeff: F) {
        self.0[wire.index] = coeff;
    }

    pub fn eval(&self, witness: &[F]) -> F {
        let mut result = F::ZERO;
        for (i, coeff) in self.0.iter().enumerate() {
            if *coeff != F::ZERO {
                result += *coeff * witness[i];
            }
        }
        result
    }

    pub fn nonzero_entries(&self) -> Vec<(usize, F)> {
        let mut result = Vec::new();
        for (i, coeff) in self.0.iter().enumerate() {
            if *coeff != F::ZERO {
                result.push((i, *coeff));
            }
        }
        result
    }
}

#[derive(Debug, Clone)]
pub struct Constraint<F: FieldExt> {
    pub A: LinearCombination<F>,
    pub B: LinearCombination<F>,
    pub C: LinearCombination<F>,
}

impl<F: FieldExt> Constraint<F> {
    pub fn new(max_terms: usize) -> Self {
        Constraint {
            A: LinearCombination::new(max_terms),
            B: LinearCombination::new(max_terms),
            C: LinearCombination::new(max_terms),
        }
    }

    pub fn is_sat(&self, witness: &[F]) -> bool {
        let lhs = self.A.eval(witness) * self.B.eval(witness);
        let rhs = self.C.eval(witness);
        lhs == rhs
    }
}

#[derive(Clone, PartialEq)]
enum Phase {
    Uninitialized,
    CounterWires,
    Synthesize,
}

#[derive(Clone, PartialEq)]
enum Mode {
    Unselected,
    WitnessGen,
    ConstraintsGen,
}

#[derive(Clone)]
pub struct ConstraintSystem<F: FieldExt> {
    pub wires: Vec<F>,
    pub constraints: Vec<Constraint<F>>,
    pub next_priv_wire: usize,
    pub next_pub_wire: usize,
    next_wire_id: usize,
    phase: Phase,
    mode: Mode,
    num_total_wires: Option<usize>,
    num_pub_inputs: Option<usize>,
    num_priv_inputs: Option<usize>,
    pub_wires: Vec<usize>,
}

impl<F: FieldExt> ConstraintSystem<F> {
    pub const fn new() -> Self {
        ConstraintSystem {
            wires: vec![],
            constraints: Vec::new(),
            next_priv_wire: 0,
            next_pub_wire: 0,
            next_wire_id: 1,
            phase: Phase::Uninitialized,
            mode: Mode::Unselected,
            num_total_wires: None,
            num_priv_inputs: None,
            num_pub_inputs: None,
            pub_wires: vec![],
        }
    }

    pub fn alloc_wire(&mut self) -> Wire<F> {
        let wire = if self.phase == Phase::CounterWires {
            self.num_total_wires = self.num_total_wires.map_or(Some(2), |x| Some(x + 1));
            Wire::new(self.next_wire_id, 0) // Set the index to 0 for now
        } else {
            let wire_index = if self.pub_wires.contains(&self.next_wire_id) {
                // If the next wire is a exposed later, allocate a public wire
                self.next_pub_wire += 1;
                self.next_pub_wire
            } else {
                self.next_priv_wire += 1;
                self.next_priv_wire
            };
            Wire::new(self.next_wire_id, wire_index)
        };

        self.next_wire_id += 1;
        wire
    }

    // Allocate a private value and return the index of the allocated wire
    pub fn alloc_priv_input(&mut self) -> Wire<F> {
        let wire = self.alloc_wire();
        if self.phase == Phase::CounterWires {
            self.num_priv_inputs = self.num_priv_inputs.map_or(Some(1), |x| Some(x + 1));
            wire
        } else {
            wire
        }
    }

    pub fn alloc_pub_input(&mut self) -> Wire<F> {
        let wire = if self.phase == Phase::CounterWires {
            self.num_total_wires = self.num_total_wires.map_or(Some(2), |x| Some(x + 1));
            self.num_pub_inputs = self.num_pub_inputs.map_or(Some(1), |x| Some(x + 1));
            Wire::new(self.next_wire_id, 0) // Set the index to 0 for now
        } else if self.phase == Phase::Synthesize {
            self.next_pub_wire += 1;
            Wire::new(self.next_wire_id, self.next_pub_wire)
        } else {
            panic!("Constraint system is't initialized");
        };

        self.next_wire_id += 1;
        wire
    }

    pub fn expose_public(&mut self, wire: Wire<F>) {
        if self.phase == Phase::CounterWires {
            self.num_pub_inputs = self.num_pub_inputs.map_or(Some(1), |x| Some(x + 1));
            // We do need to have a count of the wires
            // so we know which wires to expose
            self.pub_wires.push(wire.id);
        }
    }

    // The value "1" is always allocated at index 0 of the wires
    pub fn one() -> Wire<F> {
        Wire::new(0, 0)
    }

    // Return the constraint that enforces all of the additions.
    fn addition_con(&mut self) -> &mut Constraint<F> {
        if self.constraints.is_empty() {
            self.constraints.push(Constraint::new(self.z_len()));
            &mut self.constraints[0]
        } else {
            &mut self.constraints[0]
        }
    }

    pub fn add(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w3 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3.index] = self.wires[w1.index] + self.wires[w2.index];
            } else {
                let con = self.addition_con();
                con.A.increment_coeff(w1);
                con.A.increment_coeff(w2);
                con.B.set_coeff(Self::one(), F::ONE);
                con.C.increment_coeff(w3);
            }
        }

        w3
    }

    pub fn add_const(&mut self, w1: Wire<F>, c: F) -> Wire<F> {
        let w2 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w2.index] = self.wires[w1.index] + c;
            } else {
                let con = self.addition_con();
                con.A.increment_coeff(w1);
                con.A.increment_coeff_by(Self::one(), c);
                con.B.set_coeff(Self::one(), F::ONE);
                con.C.increment_coeff(w2);
            }
        }

        w2
    }

    pub fn neg(&mut self, w: Wire<F>) -> Wire<F> {
        let w2 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                self.wires[w2.index] = -self.wires[w.index];
            } else {
                let mut constraint = Constraint::new(self.z_len());
                constraint.A.increment_coeff(w);
                constraint.B.set_coeff(Self::one(), -F::ONE);
                constraint.C.increment_coeff(w2);
            }
        }

        w2
    }

    pub fn sub(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        w1.add(w2.neg(self), self)
    }

    pub fn sub_const(&mut self, w1: Wire<F>, c: F) -> Wire<F> {
        w1.add_const(-c, self)
    }

    pub fn mul(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w3 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3.index] = self.wires[w1.index] * self.wires[w2.index];
            } else {
                // Enforce w1 * w2 = w3 in current last constraints
                let mut constraint = Constraint::new(self.z_len());
                // Support enforcing coefficients greater than 1
                constraint.A.set_coeff(w1, F::ONE);
                constraint.B.set_coeff(w2, F::ONE);
                constraint.C.set_coeff(w3, F::ONE);

                self.constraints.push(constraint);
            }
        }

        w3
    }

    pub fn mul_const(&mut self, w1: Wire<F>, c: F) -> Wire<F> {
        let w3 = self.alloc_wire();
        // Enforce w1 * c = w3
        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3.index] = self.wires[w1.index] * c;
            } else {
                let mut constraint = Constraint::new(self.z_len());
                constraint.A.set_coeff(w1, c);
                constraint.B.set_coeff(Self::one(), F::ONE);
                constraint.C.set_coeff(w3, F::ONE);

                self.constraints.push(constraint);
            }
        }

        w3
    }

    pub fn square(&mut self, w: Wire<F>) -> Wire<F> {
        self.mul(w, w)
    }

    pub fn div(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w2_inv = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                let inv = self.wires[w2.index].invert();
                if inv.is_none().into() {
                    panic!("Division by zero at {} / {}", w1.label, w2.label);
                }

                self.wires[w2_inv.index] = inv.unwrap();
            }
        }

        let w3 = self.mul(w1, w2_inv);
        let w2_mul_w2_inv = self.mul(w2, w2_inv);
        self.assert_equal(w2_mul_w2_inv, Self::one());

        w3
    }

    // If w2 is zero, the output is assigned to zero.
    pub fn div_or_zero(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w2_inv = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                if self.wires[w2.index] == F::ZERO {
                    self.wires[w2_inv.index] = F::ZERO;
                } else {
                    self.wires[w2_inv.index] = self.wires[w2.index].invert().unwrap();
                }
            }
        }

        let w3 = self.mul(w1, w2_inv);
        let w2_mul_w2_inv = self.mul(w2, w2_inv);

        let conditional = w2.is_zero(self).not(self);
        self.assert_equal(w2_mul_w2_inv, conditional);

        w3
    }

    pub fn assert_equal(&mut self, w1: Wire<F>, w2: Wire<F>) {
        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                let assigned_w1 = self.wires[w1.index];
                let assigned_w2 = self.wires[w2.index];

                if assigned_w1 != assigned_w2 {
                    panic!(
                        "{:?} != {:?} \n 
                        for {} == {}
                            ",
                        assigned_w1,
                        assigned_w2,
                        w1.label(),
                        w2.label()
                    );
                }
            } else {
                let mut constraint = Constraint::new(self.z_len());

                // W1 * 1 = W2
                constraint.A.set_coeff(w1, F::ONE);
                constraint.B.set_coeff(Self::one(), F::ONE);
                constraint.C.set_coeff(w2, F::ONE);

                self.constraints.push(constraint);
            }
        }
    }

    pub fn is_equal(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        w1.sub(w2, self).is_zero(self)
    }

    pub fn assert_zero(&mut self, w: Wire<F>) {
        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                let assigned_w = self.wires[w.index];
                if assigned_w != F::ZERO {
                    panic!("{:?} should be zero but is {:?}", w.label(), assigned_w);
                }
            } else {
                let mut constraint = Constraint::new(self.z_len());

                // W * W = 0
                constraint.A.set_coeff(w, F::ONE);
                constraint.B.set_coeff(w, F::ONE);
                constraint.C.set_coeff(Self::one(), F::ZERO);

                self.constraints.push(constraint);
            }
        }
    }

    // Taking the same approach as the IsZero template form circomlib
    pub fn is_zero(&mut self, w: Wire<F>) -> Wire<F> {
        let inv = self.alloc_wire();

        if self.mode == Mode::WitnessGen {
            let assigned_w = self.wires[w.index];
            if assigned_w == F::ZERO {
                self.wires[inv.index] = F::ZERO;
            } else {
                self.wires[inv.index] = assigned_w.invert().unwrap();
            };
        }

        // out = -w * inv + 1
        let out = w.neg(self).mul(inv, self).add_const(F::ONE, self);

        // Assert out * w == 0
        out.mul(w, self).assert_zero(self);
        out
    }

    pub fn not(&mut self, w: Wire<F>) -> Wire<F> {
        if self.mode == Mode::WitnessGen {
            let assigned_w = self.wires[w.index];
            if assigned_w == F::ZERO || assigned_w == F::ONE {
                println!(
                    "Wire '{}' should be binary, but is {:?}",
                    w.label(),
                    assigned_w
                );
            }
        }

        w.neg(self).add_const(F::ONE, self)
    }

    // Return the number of private wires
    pub fn num_vars(&self) -> usize {
        if self.num_total_wires.is_none() {
            panic!("Number of wires not yet counted");
        }

        self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap_or(0) - 1
    }

    pub fn priv_wires_offset(&self) -> usize {
        if self.num_total_wires.is_none() {
            panic!("Number of wires not yet counted");
        }

        let num_priv_wires = self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap_or(0) - 1;
        num_priv_wires.next_power_of_two()
    }

    pub fn z_len(&self) -> usize {
        self.priv_wires_offset() * 2
    }

    // Generate the witness
    pub fn gen_witness<S: Fn(&mut ConstraintSystem<F>)>(
        &mut self,
        synthesizer: S,
        pub_inputs: &[F],
        priv_inputs: &[F],
    ) -> Vec<F> {
        if self.num_total_wires.is_none() {
            // Count the number of wires only if it hasn't been done yet
            self.phase = Phase::CounterWires;
            (synthesizer)(self);
        }

        self.wires
            .extend_from_slice(&[&[F::ONE], pub_inputs].concat());
        self.wires.resize(self.priv_wires_offset(), F::ZERO);
        self.wires.extend_from_slice(priv_inputs);
        self.wires.resize(self.z_len(), F::ZERO);

        self.start_synthesize();
        self.mode = Mode::WitnessGen;
        (synthesizer)(self);
        self.end_synthesize();

        self.wires[self.priv_wires_offset()..].to_vec()
    }

    // Produce an instance of the `R1CS` struct
    pub fn gen_constraints<S: Fn(&mut ConstraintSystem<F>)>(&mut self, synthesizer: S) -> R1CS<F> {
        if self.num_total_wires.is_none() {
            // Count the number of wires only if it hasn't been done yet
            self.phase = Phase::CounterWires;
            (synthesizer)(self);
        }

        self.start_synthesize();
        self.mode = Mode::ConstraintsGen;
        (synthesizer)(self);
        self.end_synthesize();

        let mut A_entries = vec![];
        let mut B_entries = vec![];
        let mut C_entries = vec![];

        for (i, constraint) in self.constraints.iter().enumerate() {
            for (j, coeff) in &constraint.A.nonzero_entries() {
                A_entries.push(SparseMatrixEntry {
                    row: i,
                    col: *j,
                    val: *coeff,
                });
            }

            for (j, coeff) in &constraint.B.nonzero_entries() {
                B_entries.push(SparseMatrixEntry {
                    row: i,
                    col: *j,
                    val: *coeff,
                });
            }

            for (j, coeff) in &constraint.C.nonzero_entries() {
                C_entries.push(SparseMatrixEntry {
                    row: i,
                    col: *j,
                    val: *coeff,
                });
            }
        }

        let num_cols = self.z_len();

        let A = Matrix {
            entries: A_entries,
            num_rows: self.constraints.len(),
            num_cols,
        };

        let B = Matrix {
            entries: B_entries,
            num_rows: self.constraints.len(),
            num_cols,
        };

        let C = Matrix {
            entries: C_entries,
            num_rows: self.constraints.len(),
            num_cols,
        };

        R1CS {
            A,
            B,
            C,
            num_input: self.num_pub_inputs.unwrap_or(0),
            num_vars: self.num_vars(),
        }
    }

    fn start_synthesize(&mut self) {
        self.next_wire_id = 1;
        self.next_priv_wire = self.priv_wires_offset() - 1;
        self.phase = Phase::Synthesize;
    }

    fn end_synthesize(&mut self) {
        self.phase = Phase::Uninitialized;
        self.next_wire_id = 1;
        self.next_priv_wire = 0;
        self.next_pub_wire = 0;
    }

    pub fn is_sat(
        &mut self,
        witness: &[F],
        synthesizer: impl Fn(&mut ConstraintSystem<F>),
    ) -> bool {
        self.gen_constraints(synthesizer);

        for (i, constraint) in self.constraints.iter().enumerate() {
            if !constraint.is_sat(witness) {
                println!("Constraint {} not satisfied", i);
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::mock_circuit;

    type F = shockwave_plus::halo2curves::secp256k1::Fp;

    #[test]
    fn test_phase_count_wires() {
        let (synthesizer, _, _, _) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        cs.phase = Phase::CounterWires;
        (synthesizer)(&mut cs);

        assert_eq!(cs.num_total_wires, Some(8));
        assert_eq!(cs.num_priv_inputs, Some(1));
        assert_eq!(cs.num_pub_inputs, Some(3));
    }

    #[test]
    fn test_gen_witness() {
        let (synthesizer, pub_inputs, priv_inputs, expected_witness) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let witness = cs.gen_witness(synthesizer, &pub_inputs, &priv_inputs);

        assert_eq!(cs.priv_wires_offset(), 4);
        assert_eq!(witness, expected_witness);
    }

    #[test]
    fn test_gen_constraints() {
        let (synthesizer, pub_inputs, priv_inputs, expected_witness) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.gen_constraints(&synthesizer);

        // The number of columns should equal the number of wires
        assert_eq!(cs.z_len() / 2, expected_witness.len());
        assert_eq!(r1cs.A.num_cols, cs.z_len());
        assert_eq!(r1cs.B.num_cols, cs.z_len());
        assert_eq!(r1cs.C.num_cols, cs.z_len());
    }

    #[test]
    fn test_satisfiability() {
        let (synthesizer, pub_inputs, priv_inputs, expected_witness) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.gen_constraints(&synthesizer);
        let witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        let z = R1CS::construct_z(&witness, &pub_inputs);
        assert!(cs.is_sat(&z, &synthesizer));
        assert!(r1cs.is_sat(&witness, &pub_inputs));
    }
}
