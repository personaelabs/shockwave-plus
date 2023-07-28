use crate::FieldExt;

pub fn gen_config<F: FieldExt>(num_cols: usize) -> Vec<Vec<F>> {
    assert!(num_cols.is_power_of_two());
    let expansion_factor = 2;
    let codeword_len = num_cols * expansion_factor;
    let domain = (0..codeword_len)
        .map(|i| F::from((i + 3) as u64))
        .collect::<Vec<F>>();

    let mut domain_powers = Vec::with_capacity(codeword_len);
    for eval_at in domain {
        let mut powers_i = vec![F::ONE];
        for j in 0..(num_cols - 1) {
            powers_i.push(powers_i[j] * eval_at);
        }
        domain_powers.push(powers_i);
    }

    domain_powers
}
