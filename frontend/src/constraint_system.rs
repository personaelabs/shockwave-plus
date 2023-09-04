use ark_std::{end_timer, start_timer};
use shockwave_plus::ark_ff::PrimeField;
use shockwave_plus::{Matrix, SparseMatrixEntry, R1CS};
use std::cmp::max;
use std::collections::BTreeMap;

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

    pub fn else_then(&self, out: Wire<F>) -> Wire<F> {
        self.undecided.mul(out).add(self.out)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Wire<F: PrimeField> {
    id: usize,
    index: usize,
    label: &'static str,
    cs: *mut ConstraintSystem<F>,
}

impl<F: PrimeField> Wire<F> {
    pub fn new(id: usize, index: usize, cs: *mut ConstraintSystem<F>) -> Self {
        Wire {
            id,
            index,
            label: "",
            cs,
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
    A_first: BTreeMap<usize, F>,
    B_first: BTreeMap<usize, F>,
    C_first: BTreeMap<usize, F>,
    A: BTreeMap<usize, F>,
    B: BTreeMap<usize, F>,
    C: BTreeMap<usize, F>,
    A_nonzero_coeffs: Vec<Vec<usize>>,
    B_nonzero_coeffs: Vec<Vec<usize>>,
    C_nonzero_coeffs: Vec<Vec<usize>>,
    pub next_priv_wire: usize,
    pub next_pub_wire: usize,
    next_constraint: usize,
    next_wire_id: usize,
    phase: Phase,
    mode: Mode,
    num_constraints: Option<usize>,
    num_total_wires: Option<usize>,
    num_pub_inputs: Option<usize>,
    num_priv_inputs: Option<usize>,
    pub_wires: Vec<usize>,
    z_len: usize,
}

impl<F: PrimeField> ConstraintSystem<F> {
    pub const fn new() -> Self {
        ConstraintSystem {
            wires: vec![],
            A_first: BTreeMap::new(),
            B_first: BTreeMap::new(),
            C_first: BTreeMap::new(),
            A: BTreeMap::new(),
            B: BTreeMap::new(),
            C: BTreeMap::new(),
            A_nonzero_coeffs: Vec::new(),
            B_nonzero_coeffs: Vec::new(),
            C_nonzero_coeffs: Vec::new(),
            next_priv_wire: 0,
            next_pub_wire: 0,
            next_wire_id: 1,
            phase: Phase::Uninitialized,
            mode: Mode::Unselected,
            num_total_wires: None,
            num_priv_inputs: None,
            num_pub_inputs: None,
            num_constraints: None,
            pub_wires: vec![],
            next_constraint: 1,
            z_len: 0,
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

    const fn one_wire_index() -> usize {
        0
    }

    fn next_constraint_offset(&mut self) -> usize {
        let next_constraint = self.next_constraint;
        self.next_constraint += 1;

        next_constraint * self.z_len
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

    fn increment_tree_val(tree: &mut BTreeMap<usize, F>, key: usize, val: F) {
        if let Some(v) = tree.get(&key) {
            tree.insert(key, *v + val);
        } else {
            tree.insert(key, F::ONE);
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
                Self::increment_tree_val(&mut self.A_first, w1.index, F::ONE);
                Self::increment_tree_val(&mut self.A_first, w2.index, F::ONE);
                Self::increment_tree_val(&mut self.C_first, w3.index, F::ONE);
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
                Self::increment_tree_val(&mut self.A_first, w1.index, F::ONE);
                Self::increment_tree_val(&mut self.A_first, Self::one_wire_index(), c);
                Self::increment_tree_val(&mut self.C_first, w2.index, F::ONE);
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

                let con = self.next_constraint_offset();

                let a_key = con + w.index;
                let b_key = con + Self::one_wire_index();
                let c_key = con + w2.index;

                self.A.insert(a_key, F::ONE);
                self.B.insert(b_key, -F::ONE);
                self.C.insert(c_key, F::ONE);

                self.A_nonzero_coeffs.push(vec![w.index]);
                self.B_nonzero_coeffs.push(vec![Self::one_wire_index()]);
                self.C_nonzero_coeffs.push(vec![w2.index]);
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
                let con = self.next_constraint_offset();

                let a_key = con + w1.index;
                let b_key = con + w2.index;
                let c_key = con + w3.index;

                self.A.insert(a_key, F::ONE);
                self.B.insert(b_key, F::ONE);
                self.C.insert(c_key, F::ONE);

                self.A_nonzero_coeffs.push(vec![w1.index]);
                self.B_nonzero_coeffs.push(vec![w2.index]);
                self.C_nonzero_coeffs.push(vec![w3.index]);
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

                // w1 * w2 - w3 = 0
                let con = self.next_constraint_offset();

                let a_key = con + w1.index;
                let b_key = con + Self::one_wire_index();
                let c_key = con + w3.index;

                self.A.insert(a_key, c);
                self.B.insert(b_key, F::ONE);
                self.C.insert(c_key, F::ONE);

                self.A_nonzero_coeffs.push(vec![w1.index]);
                self.B_nonzero_coeffs.push(vec![Self::one_wire_index()]);
                self.C_nonzero_coeffs.push(vec![w3.index]);
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
                let con = self.next_constraint_offset();

                // W1 * 1 == w2

                let a_key = con + w1.index;
                let b_key = con + Self::one_wire_index();
                let c_key = con + w2.index;

                self.A.insert(a_key, F::ONE);
                self.B.insert(b_key, F::ONE);
                self.C.insert(c_key, F::ONE);

                self.A_nonzero_coeffs.push(vec![w1.index]);
                self.B_nonzero_coeffs.push(vec![Self::one_wire_index()]);
                self.C_nonzero_coeffs.push(vec![w2.index]);
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
                // W * W = 0

                let con = self.next_constraint_offset();

                let a_key = con + w.index;
                let b_key = con + w.index;
                let c_key = con + Self::one_wire_index();

                self.A.insert(a_key, F::ONE);
                self.B.insert(b_key, F::ONE);
                self.C.insert(c_key, F::ZERO);

                self.A_nonzero_coeffs.push(vec![w.index]);
                self.B_nonzero_coeffs.push(vec![w.index]);
                self.C_nonzero_coeffs.push(vec![Self::one_wire_index()]);
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

        // Check the number of public inputs
        if self.num_pub_inputs.unwrap() != pub_inputs.len() {
            panic!(
                "Number of public inputs does not match. Expected {}, got {}",
                self.num_pub_inputs.unwrap(),
                pub_inputs.len()
            );
        }

        // Check the number of private inputs
        if self.num_priv_inputs.unwrap() != priv_inputs.len() {
            panic!(
                "Number of private inputs does not match. Expected {}, got {}",
                self.num_priv_inputs.unwrap(),
                priv_inputs.len()
            );
        }

        self.wires
            .extend_from_slice(&[&[F::ONE], pub_inputs].concat());
        self.wires.resize(self.priv_wires_offset(), F::ZERO);
        self.wires.extend_from_slice(priv_inputs);
        self.wires.resize(self.z_len(), F::ZERO);

        self.start_synthesize(Mode::WitnessGen);
        (synthesizer)(self);
        self.end_synthesize();

        self.wires[self.priv_wires_offset()..].to_vec()
    }

    pub fn meta<S: Fn(&mut ConstraintSystem<F>)>(&mut self, synthesizer: S) -> CircuitMeta {
        self.gen_constraints(synthesizer);

        CircuitMeta {
            num_pub_inputs: self.num_pub_inputs.unwrap_or(0),
            num_priv_inputs: self.num_priv_inputs.unwrap_or(0),
            num_constraints: self.num_constraints.unwrap(),
            num_variables: self.num_total_wires.unwrap() - self.num_pub_inputs.unwrap_or(0) - 1,
        }
    }

    pub fn gen_constraints<S: Fn(&mut ConstraintSystem<F>)>(&mut self, synthesizer: S) {
        let count_wires_timer = start_timer!(|| "Counting wires");
        if self.num_total_wires.is_none() {
            // Count the number of wires only if it hasn't been done yet
            self.phase = Phase::CounterWires;
            (synthesizer)(self);
        }
        end_timer!(count_wires_timer);

        let gen_constraints_timer = start_timer!(|| "Generating constraints");
        self.start_synthesize(Mode::ConstraintsGen);
        (synthesizer)(self);
        self.end_synthesize();
        end_timer!(gen_constraints_timer);
    }

    pub fn to_r1cs<S: Fn(&mut ConstraintSystem<F>)>(&mut self, synthesizer: S) -> R1CS<F> {
        self.gen_constraints(synthesizer);

        let constructing_r1cs = start_timer!(|| "Constructing R1CS");

        let mut A_entries = vec![];
        let mut B_entries = vec![];
        let mut C_entries = vec![];

        // First constraint

        for coeff_i in self.A_first.keys() {
            A_entries.push(SparseMatrixEntry {
                row: 0,
                col: *coeff_i,
                val: *self.A_first.get(coeff_i).unwrap(),
            });
        }

        for coeff_i in self.B_first.keys() {
            B_entries.push(SparseMatrixEntry {
                row: 0,
                col: *coeff_i,
                val: *self.B_first.get(coeff_i).unwrap(),
            });
        }

        for coeff_i in self.C_first.keys() {
            C_entries.push(SparseMatrixEntry {
                row: 0,
                col: *coeff_i,
                val: *self.C_first.get(coeff_i).unwrap(),
            });
        }

        for con in 1..self.num_constraints.unwrap() {
            let offset = con * self.z_len;
            for coeff_i in &self.A_nonzero_coeffs[con - 1] {
                A_entries.push(SparseMatrixEntry {
                    row: con,
                    col: *coeff_i,
                    val: *self.A.get(&(coeff_i + offset)).unwrap(),
                });
            }

            for coeff_i in &self.B_nonzero_coeffs[con - 1] {
                B_entries.push(SparseMatrixEntry {
                    row: con,
                    col: *coeff_i,
                    val: *self.B.get(&(coeff_i + offset)).unwrap(),
                });
            }

            for coeff_i in &self.C_nonzero_coeffs[con - 1] {
                C_entries.push(SparseMatrixEntry {
                    row: con,
                    col: *coeff_i,
                    val: *self.C.get(&(coeff_i + offset)).unwrap(),
                });
            }
        }

        let num_cols = self.z_len;
        let num_rows = self.num_constraints.unwrap();

        let A = Matrix {
            entries: A_entries,
            num_rows,
            num_cols,
        };

        let B = Matrix {
            entries: B_entries,
            num_rows,
            num_cols,
        };

        let C = Matrix {
            entries: C_entries,
            num_rows,
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

    fn start_synthesize(&mut self, mode: Mode) {
        self.next_wire_id = 1;
        self.next_priv_wire = self.priv_wires_offset() - 1;
        self.phase = Phase::Synthesize;

        // We assume that the number of constraints has been counted
        self.z_len = self.z_len();

        // The value `one` in the first row of the B matrix is always enabled.
        self.B_first.insert(Self::one_wire_index(), F::ONE);

        self.mode = mode;
    }

    fn end_synthesize(&mut self) {
        self.phase = Phase::Uninitialized;
        self.next_wire_id = 1;
        self.next_priv_wire = 0;
        self.next_pub_wire = 0;

        if self.mode == Mode::ConstraintsGen {
            self.num_constraints = Some(self.next_constraint)
        }

        self.mode = Mode::Unselected;
    }

    pub fn is_sat<S: Fn(&mut ConstraintSystem<F>)>(
        &mut self,
        witness: &[F],
        public_input: &[F],
        synthesizer: S,
    ) -> bool {
        let z = R1CS::construct_z(witness, public_input);

        self.gen_constraints(synthesizer);

        // Check the first constraint, which encodes all the additions

        let A_first_eval = self
            .A_first
            .keys()
            .map(|coeff| z[*coeff] * self.A_first.get(coeff).unwrap())
            .sum::<F>();

        let B_first_eval = self
            .B_first
            .keys()
            .map(|coeff| z[*coeff] * self.B_first.get(coeff).unwrap())
            .sum::<F>();

        let C_first_eval = self
            .C_first
            .keys()
            .map(|coeff| z[*coeff] * self.C_first.get(coeff).unwrap())
            .sum::<F>();

        if A_first_eval * B_first_eval != C_first_eval {
            println!("First constraint not satisfied");
            return false;
        }

        // Check rest of the constraints

        for con in 1..self.num_constraints.unwrap() {
            let offset = con * self.z_len;
            let A_eval: F = self.A_nonzero_coeffs[con - 1]
                .iter()
                .map(|coeff| z[*coeff] * self.A.get(&(coeff + offset)).unwrap())
                .sum();

            let B_eval: F = self.B_nonzero_coeffs[con - 1]
                .iter()
                .map(|coeff| z[*coeff] * self.B.get(&(coeff + offset)).unwrap())
                .sum();

            let C_eval: F = self.C_nonzero_coeffs[con - 1]
                .iter()
                .map(|coeff| z[*coeff] * self.C.get(&(coeff + offset)).unwrap())
                .sum();

            if A_eval * B_eval != C_eval {
                println!("Constraint {} not satisfied", con);
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
    fn test_to_r1cs() {
        let (synthesizer, _, _, expected_witness) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.to_r1cs(&synthesizer);

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

        let r1cs = cs.to_r1cs(&synthesizer);
        let witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        assert!(cs.is_sat(&witness, &pub_inputs, &synthesizer));
        assert!(r1cs.is_sat(&witness, &pub_inputs));
    }

    #[test]
    fn test_invalid_witness() {
        let (synthesizer, pub_inputs, priv_inputs, _) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.to_r1cs(&synthesizer);
        let mut witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        witness[0] += F::from(1u32);

        assert_eq!(cs.is_sat(&witness, &pub_inputs, &synthesizer), false);
        assert_eq!(r1cs.is_sat(&witness, &pub_inputs), false);
    }

    #[test]
    fn test_invalid_pub_input() {
        let (synthesizer, mut pub_inputs, priv_inputs, _) = mock_circuit();
        let mut cs = ConstraintSystem::<F>::new();

        let r1cs = cs.to_r1cs(&synthesizer);
        let witness = cs.gen_witness(&synthesizer, &pub_inputs, &priv_inputs);

        pub_inputs[0] += F::from(1u32);

        assert_eq!(cs.is_sat(&witness, &pub_inputs, &synthesizer), false);
        assert_eq!(r1cs.is_sat(&witness, &pub_inputs), false);
    }
}
