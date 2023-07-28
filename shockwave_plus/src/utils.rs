use crate::FieldExt;

// Returns a vector of vectors of length m, where each vector is a boolean vector (little endian)
pub fn boolean_hypercube<F: FieldExt>(m: usize) -> Vec<Vec<F>> {
    let n = 2usize.pow(m as u32);

    let mut boolean_hypercube = Vec::<Vec<F>>::with_capacity(n);

    for i in 0..n {
        let mut tmp = Vec::with_capacity(m);
        for j in 0..m {
            let i_b = F::from((i >> j & 1) as u64);
            tmp.push(i_b);
        }
        boolean_hypercube.push(tmp);
    }

    boolean_hypercube
}
