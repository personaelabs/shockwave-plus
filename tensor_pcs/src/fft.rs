use crate::FieldExt;

pub fn fft<F>(coeffs: &[F], domain: &[F]) -> Vec<F>
where
    F: FieldExt,
{
    debug_assert_eq!(coeffs.len(), domain.len());
    if coeffs.len() == 1 {
        return coeffs.to_vec();
    }

    // Split into evens and odds
    let L = coeffs
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 2 == 0)
        .map(|(_, x)| *x)
        .collect::<Vec<F>>();

    let R = coeffs
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 2 == 1)
        .map(|(_, x)| *x)
        .collect::<Vec<F>>();

    // Square the domain values
    let domain_squared: Vec<F> = (0..(domain.len() / 2)).map(|i| domain[i * 2]).collect();

    let fft_e = fft(&L, &domain_squared);
    let fft_o = fft(&R, &domain_squared);

    let mut evals_L = Vec::with_capacity(coeffs.len() / 2);
    let mut evals_R = Vec::with_capacity(coeffs.len() / 2);
    for i in 0..(coeffs.len() / 2) {
        // We can use the previous evaluations to create a list of evaluations
        // of the domain
        evals_L.push(fft_e[i] + fft_o[i] * domain[i]);
        evals_R.push(fft_e[i] - fft_o[i] * domain[i]);
    }

    evals_L.extend(evals_R);
    return evals_L;
}

#[cfg(test)]
mod tests {
    use halo2curves::ff::Field;
    use halo2curves::ff::PrimeField;
    use halo2curves::pasta::Fp;

    // Test the fft function by running the inverse fft
    fn ifft<F: FieldExt + Field>(domain: &[F], evals: &[F]) -> Vec<F> {
        let mut coeffs = vec![];
        let len_mod_inv = F::from(domain.len() as u64).invert().unwrap();
        let vals = fft(&evals, &domain);

        coeffs.push(vals[0] * len_mod_inv);
        for val in vals[1..].iter().rev() {
            coeffs.push(*val * len_mod_inv);
        }

        coeffs
    }

    use super::*;
    #[test]
    fn test_fft_ifft() {
        // f(x) = 1 + 2x + 3x^2 + 4x^3
        let mut coeffs = vec![
            Fp::from(1),
            Fp::from(2),
            Fp::from(3),
            Fp::from(4),
            Fp::from(5),
            Fp::from(6),
            Fp::from(7),
            Fp::from(81),
        ];

        let mut domain = vec![];

        let root_of_unity = Fp::ROOT_OF_UNITY;

        let subgroup_order = (coeffs.len() * 2).next_power_of_two();

        coeffs.resize(subgroup_order, Fp::ZERO);

        // Generator for the subgroup with order _subgroup_order_ in the field
        let generator = root_of_unity.pow(&[
            2u32.pow(32 - ((subgroup_order as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);

        for i in 0..(subgroup_order) {
            domain.push(generator.pow(&[i as u64, 0, 0, 0]));
        }

        let mut expected_evals = vec![];

        for w in &domain {
            let mut eval = Fp::ZERO;
            for (i, coeff) in (&coeffs).iter().enumerate() {
                eval += *coeff * w.pow(&[i as u64, 0, 0, 0]);
            }
            expected_evals.push(eval);
        }

        let evals = fft(&coeffs, &domain);
        debug_assert!(evals == expected_evals);

        let recovered_coeffs = ifft(&domain, &evals);
        debug_assert!(recovered_coeffs == coeffs);
    }
}
