use shockwave_plus::{Matrix, SparseMatrixEntry, R1CS};

use crate::FieldExt;

// Stores a linear combination a_0 * w_0 + a_1 * w_1 + ... + a_{n-1} * w_{n-1}
// in a sparse vector of tuples.
#[derive(Debug, Clone)]
pub struct LinearCombination<F>(Vec<F>);

impl<F: FieldExt> LinearCombination<F> {
    pub fn new(max_terms: usize) -> Self {
        let terms = vec![F::ZERO; max_terms];
        Self(terms)
    }

    pub fn increment_coeff(&mut self, wire: usize) {
        self.0[wire] += F::ONE;
    }

    pub fn set_coeff(&mut self, wire: usize, coeff: F) {
        self.0[wire] = coeff;
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
    phase: Phase,
    mode: Mode,
    wire_labels: Vec<&'static str>,
    num_total_wires: Option<usize>,
    num_pub_inputs: Option<usize>,
    num_priv_inputs: Option<usize>,
    pub_wires: Vec<&'static str>,
}

impl<F: FieldExt> ConstraintSystem<F> {
    pub const fn new() -> Self {
        ConstraintSystem {
            wires: vec![],
            constraints: Vec::new(),
            next_priv_wire: 0,
            next_pub_wire: 0,
            phase: Phase::Uninitialized,
            mode: Mode::Unselected,
            wire_labels: vec![],
            num_total_wires: None,
            num_priv_inputs: None,
            num_pub_inputs: None,
            pub_wires: vec![],
        }
    }

    pub fn alloc_wire(&mut self, label: &'static str) -> usize {
        if self.phase == Phase::CounterWires {
            self.num_total_wires = self.num_total_wires.map_or(Some(2), |x| Some(x + 1));
            0 // Return the dummy wire
        } else {
            let wire = if self.pub_wires.contains(&label) {
                // If the next wire is a public wire, allocate a public wire
                self.next_pub_wire += 1;
                println!(
                    "Allocating pub wire {} with label {}",
                    self.next_pub_wire, label
                );
                self.next_pub_wire
            } else {
                self.next_priv_wire += 1;
                self.next_priv_wire
            };
            self.wire_labels[wire] = label;
            wire
        }
    }

    // Allocate a private value and return the index of the allocated wire
    pub fn alloc_priv_input(&mut self, label: &'static str) -> usize {
        let wire = self.alloc_wire(label);
        if self.phase == Phase::CounterWires {
            self.num_priv_inputs = self.num_priv_inputs.map_or(Some(1), |x| Some(x + 1));
            0 // Return a dummy wire index
        } else {
            wire
        }
    }

    pub fn alloc_pub_input(&mut self, label: &'static str) -> usize {
        if self.phase == Phase::CounterWires {
            // Call `alloc_wire` just to increment the counter
            self.alloc_wire(label);

            self.num_pub_inputs = self.num_pub_inputs.map_or(Some(1), |x| Some(x + 1));
            0 // Return a dummy wire index
        } else if self.phase == Phase::Synthesize {
            self.next_pub_wire += 1;
            self.next_pub_wire
        } else {
            panic!("Constraint system is't initialized");
        }
    }

    pub fn expose_public(&mut self, wire_label: &'static str) {
        if self.phase == Phase::CounterWires {
            self.num_pub_inputs = self.num_pub_inputs.map_or(Some(1), |x| Some(x + 1));
            // We do need to have a count of the wires
            // so we know which wires to expose
            self.pub_wires.push(wire_label);
        }
    }

    // The value "1" is always allocated at index 0 of the wires
    pub fn one_index() -> usize {
        0
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

    // Enforce that w1 + w2 = w3 where w1, w2, w3 are the indices of the wires
    // compiling the witness generator
    pub fn add(&mut self, w1: usize, w2: usize, label: &'static str) -> usize {
        let w3 = self.alloc_wire(label);

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3] = self.wires[w1] + self.wires[w2];
            } else {
                let con = self.addition_con();
                con.A.increment_coeff(w1);
                con.A.increment_coeff(w2);
                con.B.set_coeff(Self::one_index(), F::ONE);
                con.C.increment_coeff(w3);
            }
        }

        w3
    }

    pub fn mul(&mut self, w1: usize, w2: usize, label: &'static str) -> usize {
        let w3 = self.alloc_wire(label);

        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3] = self.wires[w1] * self.wires[w2];
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

    pub fn mul_const(&mut self, w1: usize, c: F, label: &'static str) -> usize {
        let w3 = self.alloc_wire(label);
        // Enforce w1 * c = w3
        if self.phase == Phase::Synthesize {
            if self.mode == Mode::WitnessGen {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3] = self.wires[w1] * c;
            } else {
                let mut constraint = Constraint::new(self.z_len());
                constraint.A.set_coeff(w1, c);
                constraint.B.set_coeff(Self::one_index(), F::ONE);
                constraint.C.set_coeff(w3, F::ONE);

                self.constraints.push(constraint);
            }
        }

        w3
    }

    // Return the number of private wires
    pub fn num_vars(&self) -> usize {
        if self.num_total_wires.is_none() {
            panic!("Number of wires not yet counted");
        }

        self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap() - 1
    }

    pub fn priv_wires_offset(&self) -> usize {
        if self.num_total_wires.is_none() {
            panic!("Number of wires not yet counted");
        }

        let num_priv_wires = self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap() - 1;
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
            num_input: self.num_pub_inputs.unwrap(),
            num_vars: self.num_vars(),
        }
    }

    fn start_synthesize(&mut self) {
        self.next_priv_wire = self.priv_wires_offset() - 1;
        self.wire_labels.resize(self.z_len(), "");
        self.phase = Phase::Synthesize;
    }

    fn end_synthesize(&mut self) {
        self.phase = Phase::Uninitialized;
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

    type F = tensor_pcs::halo2curves::secp256k1::Fp;

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
