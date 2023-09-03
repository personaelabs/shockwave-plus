use std::cmp::max;
use std::marker::PhantomData;

use ark_std::{end_timer, start_timer};
use shockwave_plus::{Matrix, SparseMatrixEntry, R1CS};

use shockwave_plus::ark_ff::PrimeField;

pub struct Conditional<F: PrimeField> {
    undecided: Wire<F>,
    out: Wire<F>,
}

impl<F: PrimeField> Conditional<F> {
    pub fn if_then(sel: Wire<F>, out: Wire<F>, cs: &mut ConstraintSystem<F>) -> Self {
        cs.assert_binary(sel);

        let out = sel.mul(out);

        Self {
            undecided: !sel,
            out,
        }
    }

    pub fn elif(&self, sel: Wire<F>, out: Wire<F>, cs: &mut ConstraintSystem<F>) -> Self {
        cs.assert_binary(sel);

        let this_cond = sel.mul(out);
        let out = self.undecided.mul(this_cond).add(self.out);
        let undecided = !sel & self.undecided;

        Self { undecided, out }
    }

    pub fn else_then(&self, out: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        self.undecided.mul(out).add(self.out)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Wire<F: PrimeField> {
    id: usize,
    index: usize,
    label: &'static str,
    cs: *mut ConstraintSystem<F>,
    _marker: PhantomData<F>,
}

impl<F: PrimeField> Wire<F> {
    pub fn new(id: usize, index: usize, cs: *mut ConstraintSystem<F>) -> Self {
        Wire {
            id,
            index,
            label: "",
            cs,
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

    pub fn cs(&self) -> &mut ConstraintSystem<F> {
        unsafe { &mut *self.cs as &mut ConstraintSystem<F> }
    }

    pub fn div_or_zero(&self, w: Wire<F>) -> Wire<F> {
        w.cs().div_or_zero(*self, w)
    }

    pub fn is_zero(&self) -> Wire<F> {
        self.cs().is_zero(*self)
    }

    pub fn is_equal(&self, w: Wire<F>) -> Wire<F> {
        self.cs().is_equal(*self, w)
    }

    pub fn assert_zero(&self, cs: &mut ConstraintSystem<F>) {
        cs.assert_zero(*self)
    }

    pub fn assert_equal(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) {
        cs.assert_equal(*self, w)
    }

    pub fn and(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.and(*self, w)
    }

    pub fn or(&self, w: Wire<F>, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.or(*self, w)
    }

    pub fn not(&self, cs: &mut ConstraintSystem<F>) -> Wire<F> {
        cs.not(*self)
    }

    pub fn print_val(&self) {
        if self.cs().mode == Mode::WitnessGen {
            println!(
                "{} = {:?}",
                self.label(),
                self.cs().wires[self.index].to_string()
            );
        }
    }

    pub fn val(&self, cs: &mut ConstraintSystem<F>) -> Option<F> {
        if cs.mode == Mode::WitnessGen {
            Some(cs.wires[self.index])
        } else {
            None
        }
    }
}

use std::ops::{Add, AddAssign, BitAnd, BitOr, Div, Mul, Neg, Not, Sub, SubAssign};

impl<F: PrimeField> Add<Wire<F>> for Wire<F> {
    type Output = Wire<F>;

    fn add(self, rhs: Wire<F>) -> Self::Output {
        self.cs().add(self, rhs)
    }
}

impl<F: PrimeField> AddAssign<Wire<F>> for Wire<F> {
    fn add_assign(&mut self, rhs: Wire<F>) {
        *self = self.cs().add(*self, rhs);
    }
}

impl<F: PrimeField> Sub<Wire<F>> for Wire<F> {
    type Output = Wire<F>;

    fn sub(self, rhs: Wire<F>) -> Self::Output {
        self.cs().sub(self, rhs)
    }
}

impl<F: PrimeField> SubAssign<Wire<F>> for Wire<F> {
    fn sub_assign(&mut self, rhs: Wire<F>) {
        *self = self.cs().sub(*self, rhs);
    }
}

impl<F: PrimeField> Mul<Wire<F>> for Wire<F> {
    type Output = Wire<F>;

    fn mul(self, rhs: Wire<F>) -> Self::Output {
        self.cs().mul(self, rhs)
    }
}

impl<F: PrimeField> Div<Wire<F>> for Wire<F> {
    type Output = Wire<F>;

    fn div(self, rhs: Wire<F>) -> Self::Output {
        self.cs().div(self, rhs)
    }
}

impl<F: PrimeField> Neg for Wire<F> {
    type Output = Wire<F>;

    fn neg(self) -> Self::Output {
        self.cs().neg(self)
    }
}

impl<F: PrimeField> BitAnd for Wire<F> {
    type Output = Wire<F>;

    fn bitand(self, rhs: Wire<F>) -> Self::Output {
        self.cs().and(self, rhs)
    }
}

impl<F: PrimeField> BitOr for Wire<F> {
    type Output = Wire<F>;

    fn bitor(self, rhs: Wire<F>) -> Self::Output {
        self.cs().or(self, rhs)
    }
}

impl<F: PrimeField> Not for Wire<F> {
    type Output = Wire<F>;

    fn not(self) -> Self::Output {
        self.cs().not(self)
    }
}

// Stores a linear combination a_0 * w_0 + a_1 * w_1 + ... + a_{n-1} * w_{n-1}
// in a sparse vector of tuples.
#[derive(Debug, Clone)]
pub struct LinearCombination<F>(Vec<F>);

impl<F: PrimeField> LinearCombination<F> {
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
pub struct Constraint<F: PrimeField> {
    pub A: LinearCombination<F>,
    pub B: LinearCombination<F>,
    pub C: LinearCombination<F>,
}

impl<F: PrimeField> Constraint<F> {
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

#[derive(Debug, Clone, Copy)]
pub struct CircuitMeta {
    pub num_pub_inputs: usize,
    pub num_priv_inputs: usize,
    pub num_constraints: usize,
    pub num_variables: usize,
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

// A constraint system (R1CS) that handles the followings:
// - Allocating and constraining wires (by exposing methods `add`, `mul`, etc.)
// - Generating the witness (`gen_witness`)
// - Generating the R1CS instance (`gen_constraints`)
// The circuit writer should define a "synthesizer" that takes a mutable reference
// to a `ConstraintSystem` and calls its method to allocate and constrain wires.
#[derive(Clone)]
pub struct ConstraintSystem<F: PrimeField> {
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

impl<F: PrimeField> ConstraintSystem<F> {
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

    pub fn is_witness_gen(&self) -> bool {
        self.mode == Mode::WitnessGen
    }

    fn alloc_wire(&mut self) -> Wire<F> {
        let wire = if self.phase == Phase::CounterWires {
            self.num_total_wires = self.num_total_wires.map_or(Some(2), |x| Some(x + 1));
            Wire::new(self.next_wire_id, 0, self) // Set the index to 0 for now
        } else {
            let wire_index = if self.pub_wires.contains(&self.next_wire_id) {
                // If the next wire is a exposed later, allocate a public wire
                self.next_pub_wire += 1;
                self.next_pub_wire
            } else {
                self.next_priv_wire += 1;
                self.next_priv_wire
            };
            Wire::new(self.next_wire_id, wire_index, self)
        };

        self.next_wire_id += 1;
        wire
    }

    // Allocate an unconstrained variable.
    // Use `alloc_const` to allocate a constant value.
    pub fn alloc_var(&mut self, val: F) -> Wire<F> {
        let wire = self.alloc_wire();
        if self.is_witness_gen() {
            self.wires[wire.index] = val;
        }
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

    pub fn alloc_priv_inputs(&mut self, n: usize) -> Vec<Wire<F>> {
        (0..n).map(|_| self.alloc_priv_input()).collect()
    }

    // Allocate a public input wire.
    pub fn alloc_pub_input(&mut self) -> Wire<F> {
        let wire = if self.phase == Phase::CounterWires {
            self.num_total_wires = self.num_total_wires.map_or(Some(2), |x| Some(x + 1));
            self.num_pub_inputs = self.num_pub_inputs.map_or(Some(1), |x| Some(x + 1));
            Wire::new(self.next_wire_id, 0, self) // Set the index to 0 for now
        } else if self.phase == Phase::Synthesize {
            self.next_pub_wire += 1;
            Wire::new(self.next_wire_id, self.next_pub_wire, self)
        } else {
            panic!("Constraint system is't initialized");
        };

        self.next_wire_id += 1;
        wire
    }

    // Expose a wire as a public input.
    pub fn expose_public(&mut self, wire: Wire<F>) {
        if self.phase == Phase::CounterWires {
            self.num_pub_inputs = self.num_pub_inputs.map_or(Some(1), |x| Some(x + 1));
            // We need to count wires so we know which wires to expose.
            self.pub_wires.push(wire.id);
        }
    }

    // Allocate a constant value.
    pub fn alloc_const(&mut self, c: F) -> Wire<F> {
        let one = self.one();
        self.mul_const(one, c)
    }

    // The value "1" is a
    pub fn one(&mut self) -> Wire<F> {
        Wire::new(0, 0, self)
    }

    // Return the constraint that enforces all of the additions and subtractions.
    fn addition_con(&mut self) -> &mut Constraint<F> {
        if self.constraints.is_empty() {
            self.constraints.push(Constraint::new(self.z_len()));
            &mut self.constraints[0]
        } else {
            &mut self.constraints[0]
        }
    }

    // Assert that the given wire is binary at witness generation.
    // It does NOT constraint the wire to be binary.
    fn assert_binary(&self, w: Wire<F>) {
        if self.is_witness_gen() {
            let assigned_w = self.wires[w.index];
            if assigned_w != F::ZERO && assigned_w != F::ONE {
                println!(
                    "Wire '{}' should be binary, but is {:?}",
                    w.label(),
                    assigned_w
                );
            }
        }
    }

    pub fn add(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w3 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                // If this is a witness generation call,
                // we assign the output of the gate here.
                self.wires[w3.index] = self.wires[w1.index] + self.wires[w2.index];
            } else {
                // (w1 + w2) * 1 - w3 = 0
                let one = self.one();
                let con = self.addition_con();
                con.A.increment_coeff(w1);
                con.A.increment_coeff(w2);
                con.B.set_coeff(one, F::ONE);
                con.C.increment_coeff(w3);
            }
        }

        w3
    }

    pub fn add_const(&mut self, w1: Wire<F>, c: F) -> Wire<F> {
        let w2 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                self.wires[w2.index] = self.wires[w1.index] + c;
            } else {
                // (w1 + c) * 1 - w2 = 0
                let one = self.one();
                let con = self.addition_con();
                con.A.increment_coeff(w1);
                con.A.increment_coeff_by(one, c);
                con.B.set_coeff(one, F::ONE);
                con.C.increment_coeff(w2);
            }
        }

        w2
    }

    pub fn neg(&mut self, w: Wire<F>) -> Wire<F> {
        let w2 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                self.wires[w2.index] = -self.wires[w.index];
            } else {
                // w * (1 * -1) - w2 = 0
                let mut constraint = Constraint::new(self.z_len());
                constraint.A.set_coeff(w, F::ONE);
                constraint.B.set_coeff(self.one(), -F::ONE);
                constraint.C.set_coeff(w2, F::ONE);
            }
        }

        w2
    }

    pub fn sub(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        w1.add(w2.neg())
    }

    // Subtract a constant value from a wire.
    pub fn sub_const(&mut self, w1: Wire<F>, c: F) -> Wire<F> {
        self.add_const(w1, -c)
    }

    pub fn mul(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w3 = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                self.wires[w3.index] = self.wires[w1.index] * self.wires[w2.index];
            } else {
                // w1 * w2 - w3 = 0
                let mut constraint = Constraint::new(self.z_len());
                constraint.A.set_coeff(w1, F::ONE);
                constraint.B.set_coeff(w2, F::ONE);
                constraint.C.set_coeff(w3, F::ONE);

                self.constraints.push(constraint);
            }
        }

        w3
    }

    // Multiply a wire by a constant value
    pub fn mul_const(&mut self, w1: Wire<F>, c: F) -> Wire<F> {
        let w3 = self.alloc_wire();
        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                self.wires[w3.index] = self.wires[w1.index] * c;
            } else {
                // w1 * c - w3 = 0
                let mut constraint = Constraint::new(self.z_len());
                constraint.A.set_coeff(w1, c);
                constraint.B.set_coeff(self.one(), F::ONE);
                constraint.C.set_coeff(w3, F::ONE);

                self.constraints.push(constraint);
            }
        }

        w3
    }

    pub fn square(&mut self, w: Wire<F>) -> Wire<F> {
        self.mul(w, w)
    }

    // This function will panic if the denominator is zero.
    // Use `div_or_zero` to handle division by zero.
    pub fn div(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w2_inv = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                let inv = self.wires[w2.index].inverse();
                if inv.is_none().into() {
                    panic!("Division by zero at {} / {}", w1.label(), w2.label());
                }

                self.wires[w2_inv.index] = inv.unwrap();
            }
        }

        let w3 = self.mul(w1, w2_inv);
        let w2_mul_w2_inv = self.mul(w2, w2_inv);
        let one = self.one();
        self.assert_equal(w2_mul_w2_inv, one);

        w3
    }

    // If the denominator is zero, the output is assigned to zero.
    pub fn div_or_zero(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        let w2_inv = self.alloc_wire();

        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                if self.wires[w2.index] == F::ZERO {
                    self.wires[w2_inv.index] = F::ZERO;
                } else {
                    self.wires[w2_inv.index] = self.wires[w2.index].inverse().unwrap();
                }
            }
        }

        let w3 = self.mul(w1, w2_inv);
        let w2_mul_w2_inv = self.mul(w2, w2_inv);

        let conditional = !(w2.is_zero());
        self.assert_equal(w2_mul_w2_inv, conditional);

        w3
    }

    pub fn assert_equal(&mut self, w1: Wire<F>, w2: Wire<F>) {
        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
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
                constraint.B.set_coeff(self.one(), F::ONE);
                constraint.C.set_coeff(w2, F::ONE);

                self.constraints.push(constraint);
            }
        }
    }

    pub fn is_equal(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        (w1 - w2).is_zero()
    }

    // Start a conditional block.
    // The conditional block must be ended by calling `else_then`.
    pub fn if_then(&mut self, sel: Wire<F>, out: Wire<F>) -> Conditional<F> {
        Conditional::if_then(sel, out, self)
    }

    pub fn assert_zero(&mut self, w: Wire<F>) {
        if self.phase == Phase::Synthesize {
            if self.is_witness_gen() {
                let assigned_w = self.wires[w.index];
                if assigned_w != F::ZERO {
                    panic!("{:?} should be zero but is {:?}", w.label(), assigned_w);
                }
            } else {
                let mut constraint = Constraint::new(self.z_len());

                // W * W = 0
                constraint.A.set_coeff(w, F::ONE);
                constraint.B.set_coeff(w, F::ONE);
                constraint.C.set_coeff(self.one(), F::ZERO);

                self.constraints.push(constraint);
            }
        }
    }

    // Return a binary wire that is 1 if the input wire is zero and 0 otherwise.
    pub fn is_zero(&mut self, w: Wire<F>) -> Wire<F> {
        // Taking the same approach as the IsZero template form circomlib

        let inv = self.alloc_wire();
        if self.is_witness_gen() {
            let assigned_w = self.wires[w.index];
            if assigned_w == F::ZERO {
                self.wires[inv.index] = F::ZERO;
            } else {
                self.wires[inv.index] = assigned_w.inverse().unwrap();
            };
        }

        let one = self.one();
        let out = w.neg().mul(inv) + one;
        self.assert_zero(out * w);
        out
    }

    // Returns `w1 AND w2`.
    // It does NOT constraint the input wires to be binary.
    pub fn and(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        self.assert_binary(w1);
        self.assert_binary(w2);
        w1.mul(w2)
    }

    // Returns `w1 OR w2`.
    // It does NOT constraint the input wires to be binary.
    pub fn or(&mut self, w1: Wire<F>, w2: Wire<F>) -> Wire<F> {
        (w1 + w2) - (w1 * w2)
    }

    // Returns `!w`.
    // It does NOT constraint the input wires to be binary.
    pub fn not(&mut self, w: Wire<F>) -> Wire<F> {
        self.assert_binary(w);
        let one = self.one();
        one - w
    }

    // Return the number of private wires
    pub fn num_vars(&self) -> usize {
        if self.num_total_wires.is_none() {
            panic!("Number of wires not yet counted");
        }

        self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap_or(0) - 1
    }

    fn det_priv_wires_offset(num_total_wires: usize, num_pub_inputs: usize) -> usize {
        // We do +1 to account for the wire that is always "1".
        let num_pub_wires = num_pub_inputs + 1;
        let num_priv_wires = num_total_wires - num_pub_wires;

        max(num_pub_wires, num_priv_wires).next_power_of_two()
    }

    fn priv_wires_offset(&self) -> usize {
        if self.num_total_wires.is_none() {
            panic!("Number of wires not yet counted");
        }

        Self::det_priv_wires_offset(
            self.num_total_wires.unwrap(),
            self.num_pub_inputs.unwrap_or(0),
        )
    }

    fn z_len(&self) -> usize {
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

    pub fn meta<S: Fn(&mut ConstraintSystem<F>)>(&mut self, synthesizer: S) -> CircuitMeta {
        self.gen_constraints(synthesizer);

        CircuitMeta {
            num_pub_inputs: self.num_pub_inputs.unwrap_or(0),
            num_priv_inputs: self.num_priv_inputs.unwrap_or(0),
            num_constraints: self.constraints.len(),
            num_variables: self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap_or(0) - 1,
        }
    }

    // Produce an instance of the `R1CS` struct
    pub fn gen_constraints<S: Fn(&mut ConstraintSystem<F>)>(&mut self, synthesizer: S) -> R1CS<F> {
        let count_wires_timer = start_timer!(|| "Counting wires");
        if self.num_total_wires.is_none() {
            // Count the number of wires only if it hasn't been done yet
            self.phase = Phase::CounterWires;
            (synthesizer)(self);
        }
        end_timer!(count_wires_timer);

        let gen_constraints_timer = start_timer!(|| "Generating constraints");
        self.start_synthesize();
        self.mode = Mode::ConstraintsGen;
        (synthesizer)(self);
        self.end_synthesize();
        end_timer!(gen_constraints_timer);

        let constructing_r1cs = start_timer!(|| "Constructing R1CS");
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

        end_timer!(constructing_r1cs);

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
        public_input: &[F],
        synthesizer: impl Fn(&mut ConstraintSystem<F>),
    ) -> bool {
        let z = R1CS::construct_z(witness, public_input);

        self.gen_constraints(synthesizer);

        for (i, constraint) in self.constraints.iter().enumerate() {
            if !constraint.is_sat(&z) {
                println!("Constraint {} not satisfied", i);
                return false;
            }
        }

        true
    }
}

#[macro_export]
macro_rules! init_constraint_system {
    ($field:ty) => {
        pub use std::sync::Mutex;

        static CONSTRAINT_SYSTEM: Mutex<ConstraintSystem<$field>> =
            Mutex::new(ConstraintSystem::new());
    };
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::test_utils::mock_circuit;

    type F = shockwave_plus::ark_secp256k1::Fq;

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
    fn test_private_wires_offset() {
        // Test that the constraint system correctly determines the offset
        // where private wires start

        // Case 1
        let num_total_wires = 3;
        let num_pub_inputs = 1;
        // num_priv_wires = 3 - (1 + 1) = 1
        let offset = ConstraintSystem::<F>::det_priv_wires_offset(num_total_wires, num_pub_inputs);
        assert_eq!(offset, 2);

        // Case 2
        let num_total_wires = 4;
        let num_pub_inputs = 1;
        // num_priv_wires = 4 - (1 + 1) = 2
        let offset = ConstraintSystem::<F>::det_priv_wires_offset(num_total_wires, num_pub_inputs);
        assert_eq!(offset, 2);

        // Case 3
        let num_total_wires = 4;
        let num_pub_inputs = 2;
        // num_priv_wires = 4 - (2 + 1) = 1
        let offset = ConstraintSystem::<F>::det_priv_wires_offset(num_total_wires, num_pub_inputs);
        assert_eq!(offset, 4);
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
        let (synthesizer, _, _, expected_witness) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.gen_constraints(&synthesizer);

        // The number of columns should equal the number of wires
        assert_eq!(cs.z_len() / 2, expected_witness.len());
        assert_eq!(r1cs.A.num_cols, cs.z_len());
        assert_eq!(r1cs.B.num_cols, cs.z_len());
        assert_eq!(r1cs.C.num_cols, cs.z_len());
    }

    #[test]
    fn test_valid_witness() {
        let (synthesizer, pub_inputs, priv_inputs, _) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.gen_constraints(&synthesizer);
        let witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        assert!(cs.is_sat(&witness, &pub_inputs, &synthesizer));
        assert!(r1cs.is_sat(&witness, &pub_inputs));
    }

    #[test]
    fn test_invalid_witness() {
        let (synthesizer, pub_inputs, priv_inputs, _) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.gen_constraints(&synthesizer);
        let mut witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        witness[0] += F::from(1u32);

        assert_eq!(cs.is_sat(&witness, &pub_inputs, &synthesizer), false);
        assert_eq!(r1cs.is_sat(&witness, &pub_inputs), false);
    }

    #[test]
    fn test_invalid_pub_input() {
        let (synthesizer, mut pub_inputs, priv_inputs, _) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.gen_constraints(&synthesizer);
        let witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        pub_inputs[0] += F::from(1u32);

        assert_eq!(cs.is_sat(&witness, &pub_inputs, &synthesizer), false);
        assert_eq!(r1cs.is_sat(&witness, &pub_inputs), false);
    }
}
