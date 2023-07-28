use crate::FieldExt;
use halo2curves::ff::Field;
use tensor_pcs::SparseMLPoly;

#[derive(Clone)]
pub struct SparseMatrixEntry<F: FieldExt> {
    pub row: usize,
    pub col: usize,
    pub val: F,
}

#[derive(Clone)]
pub struct Matrix<F: FieldExt> {
    pub entries: Vec<SparseMatrixEntry<F>>,
    pub num_cols: usize,
    pub num_rows: usize,
}

impl<F> Matrix<F>
where
    F: FieldExt,
{
    pub fn new(entries: Vec<SparseMatrixEntry<F>>, num_cols: usize, num_rows: usize) -> Self {
        assert!((num_cols * num_rows).is_power_of_two());
        Self {
            entries,
            num_cols,
            num_rows,
        }
    }

    pub fn mul_vector(&self, num_rows: usize, vec: &[F]) -> Vec<F> {
        let mut result = vec![F::ZERO; num_rows];
        let entries = &self.entries;
        for i in 0..entries.len() {
            let row = entries[i].row;
            let col = entries[i].col;
            let val = entries[i].val;
            result[row] += val * vec[col];
        }
        result
    }

    // Return a multilinear extension of the matrix
    // with num_vars * num_vars entries
    pub fn to_ml_extension(&self) -> SparseMLPoly<F> {
        let mut evals = Vec::with_capacity(self.entries.len());
        let entries = &self.entries;
        let num_cols = self.num_cols;
        for i in 0..entries.len() {
            let row = entries[i].row;
            let col = entries[i].col;
            let val = entries[i].val;
            evals.push(((row * num_cols) + col, val));
        }
        let ml_poly_num_vars = ((self.num_cols * self.num_rows) as f64).log2() as usize;
        let ml_poly = SparseMLPoly::new(evals, ml_poly_num_vars);
        ml_poly
    }

    /*
    pub fn fast_to_coeffs(&self, s: usize, x: F) -> Vec<F> {
        let mut result = F::ZERO;
        for entry in &self.0 {
            let row = entry.0;
            let col = entry.1;
            let val = entry.2;

            let index = row * 2usize.pow(s as u32) + col;
            // Get the degrees of the nonzero coefficients
            // Tensor product (1 - x_0)(1 - x_1)
            let base = index;
            let zero_bits = degree & !base;

            let mut zero_bit_degrees = vec![];
            for j in 0..s {
                if zero_bits & (1 << j) != 0 {
                    zero_bit_degrees.push(j);
                }
            }

            let mut term = val;
            for degree in zero_bit_degrees {
                term *= x.pow(&[base as u64, 0, 0, 0]) - x.pow(&[(degree + base) as u64, 0, 0, 0]);
            }
            result += term;
        }
        result
    }
     */

    /*
    pub fn fast_uni_eval(&self, s: usize, x: F) -> F {
        let degree = 2usize.pow(s as u32);

        let mut result = F::ZERO;
        for entry in &self.0 {
            let row = entry.0;
            let col = entry.1;
            let val = entry.2;

            let index = row * 2usize.pow(s as u32) + col;
            // Get the degrees of the nonzero coefficients
            // Tensor product (1 - x_0)(1 - x_1)
            let base = index;
            let zero_bits = degree & !base;

            let mut zero_bit_degrees = vec![];
            for j in 0..s {
                if zero_bits & (1 << j) != 0 {
                    zero_bit_degrees.push(j);
                }
            }

            let mut term = val;
            for degree in zero_bit_degrees {
                term *= x.pow(&[base as u64, 0, 0, 0]) - x.pow(&[(degree + base) as u64, 0, 0, 0]);
            }
            result += term;
        }
        result
    }
     */
}

#[derive(Clone)]
pub struct R1CS<F>
where
    F: FieldExt,
{
    pub A: Matrix<F>,
    pub B: Matrix<F>,
    pub C: Matrix<F>,
    pub public_input: Vec<F>,
    pub num_cons: usize,
    pub num_vars: usize,
    pub num_input: usize,
}

impl<F> R1CS<F>
where
    F: FieldExt,
{
    pub fn hadamard_prod(a: &[F], b: &[F]) -> Vec<F> {
        assert_eq!(a.len(), b.len());
        let mut result = vec![F::ZERO; a.len()];
        for i in 0..a.len() {
            result[i] = a[i] * b[i];
        }
        result
    }

    pub fn produce_synthetic_r1cs(
        num_cons: usize,
        num_vars: usize,
        num_input: usize,
    ) -> (Self, Vec<F>) {
        //        assert_eq!(num_cons, num_vars);
        let mut public_input = Vec::with_capacity(num_input);
        let mut witness = Vec::with_capacity(num_vars);

        for i in 0..num_input {
            public_input.push(F::from((i + 1) as u64));
        }

        for i in 0..num_vars {
            witness.push(F::from((i + 1) as u64));
        }

        let z: Vec<F> = vec![public_input.clone(), witness.clone()].concat();

        let mut A_entries: Vec<SparseMatrixEntry<F>> = vec![];
        let mut B_entries: Vec<SparseMatrixEntry<F>> = vec![];
        let mut C_entries: Vec<SparseMatrixEntry<F>> = vec![];

        for i in 0..num_cons {
            let A_col = i % num_vars;
            let B_col = (i + 1) % num_vars;
            let C_col = (i + 2) % num_vars;

            // For the i'th constraint,
            // add the value 1 at the (i % num_vars)th column of A, B.
            // Compute the corresponding C_column value so that A_i * B_i = C_i
            // we apply multiplication since the Hadamard product is computed for Az ・ Bz,

            // We only _enable_ a single variable in each constraint.
            A_entries.push(SparseMatrixEntry {
                row: i,
                col: A_col,
                val: F::ONE,
            });
            B_entries.push(SparseMatrixEntry {
                row: i,
                col: B_col,
                val: F::ONE,
            });
            C_entries.push(SparseMatrixEntry {
                row: i,
                col: C_col,
                val: (z[A_col] * z[B_col]) * z[C_col].invert().unwrap(),
            });
        }

        let A = Matrix::new(A_entries, num_vars, num_cons);

        let B = Matrix::new(B_entries, num_vars, num_cons);

        let C = Matrix::new(C_entries, num_vars, num_cons);

        (
            Self {
                A,
                B,
                C,
                public_input,
                num_cons,
                num_vars,
                num_input,
            },
            witness,
        )
    }

    pub fn is_sat(&self, witness: &Vec<F>, public_input: &Vec<F>) -> bool {
        let mut z = Vec::with_capacity(witness.len() + public_input.len() + 1);
        z.extend(public_input);
        z.extend(witness);

        let Az = self.A.mul_vector(self.num_cons, &z);
        let Bz = self.B.mul_vector(self.num_cons, &z);
        let Cz = self.C.mul_vector(self.num_cons, &z);

        Self::hadamard_prod(&Az, &Bz) == Cz
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::boolean_hypercube;

    use super::*;
    type F = halo2curves::secp256k1::Fp;
    use crate::polynomial::ml_poly::MlPoly;

    #[test]
    fn test_r1cs() {
        let num_cons = 2usize.pow(5);
        let num_vars = num_cons;
        let num_input = 0;

        let (r1cs, mut witness) = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        assert_eq!(witness.len(), num_vars);
        assert_eq!(r1cs.public_input.len(), num_input);

        assert!(r1cs.is_sat(&witness, &r1cs.public_input));

        // Should assert if the witness is invalid
        witness[0] = witness[0] + F::one();
        assert!(r1cs.is_sat(&r1cs.public_input, &witness) == false);
        witness[0] = witness[0] - F::one();

        /*
        // Should assert if the public input is invalid
        let mut public_input = r1cs.public_input.clone();
        public_input[0] = public_input[0] + F::one();
        assert!(r1cs.is_sat(&witness, &public_input) == false);
         */

        // Test MLE
        let s = (num_vars as f64).log2() as usize;
        let A_mle = r1cs.A.to_ml_extension();
        let B_mle = r1cs.B.to_ml_extension();
        let C_mle = r1cs.C.to_ml_extension();
        let Z_mle = MlPoly::new(witness);

        for c in &boolean_hypercube(s) {
            let mut eval_a = F::zero();
            let mut eval_b = F::zero();
            let mut eval_c = F::zero();
            for b in &boolean_hypercube(s) {
                let mut b_rev = b.clone();
                b_rev.reverse();
                let z_eval = Z_mle.eval(&b_rev);
                let mut eval_matrix = [b.as_slice(), c.as_slice()].concat();
                eval_matrix.reverse();
                eval_a += A_mle.eval(&eval_matrix) * z_eval;
                eval_b += B_mle.eval(&eval_matrix) * z_eval;
                eval_c += C_mle.eval(&eval_matrix) * z_eval;
            }
            let eval_con = eval_a * eval_b - eval_c;
            assert_eq!(eval_con, F::zero());
        }
    }

    /*
    #[test]
    fn test_fast_uni_eval() {
        let (r1cs, _) = R1CS::<F>::produce_synthetic_r1cs(8, 8, 0);

        let eval_at = F::from(33);
        let result = r1cs.A.fast_uni_eval(r1cs.num_vars, eval_at);
        println!("result: {:?}", result);
    }
     */
}
