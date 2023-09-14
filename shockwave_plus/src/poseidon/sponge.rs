use crate::{poseidon::Poseidon, FieldGC, PoseidonCurve};
use std::result::Result;
use tiny_keccak::{Hasher, Keccak};

use super::CAPACITY;

#[derive(Clone, PartialEq, Debug)]
pub enum SpongeOp {
    Absorb(usize),
    Squeeze(usize),
}

#[derive(Clone)]
pub struct IOPattern(pub Vec<SpongeOp>);

impl IOPattern {
    pub const fn new(io_pattern: Vec<SpongeOp>) -> Self {
        Self(io_pattern)
    }
}

// Implements SAFE (Sponge API for Field Elements): https://hackmd.io/bHgsH6mMStCVibM_wYvb2w
#[derive(Clone)]
pub struct PoseidonSponge<F: FieldGC, const W: usize> {
    pub absorb_pos: usize,
    pub squeeze_pos: usize,
    pub io_count: usize,
    pub io_pattern: IOPattern,
    pub rate: usize,
    pub capacity: usize,
    pub tag: F,
    poseidon: Poseidon<F, W>,
}

impl<F: FieldGC, const WIDTH: usize> PoseidonSponge<F, WIDTH> {
    pub fn new(domain_separator: &[u8], curve: PoseidonCurve, io_pattern: IOPattern) -> Self {
        // Parse the constants from string

        let tag = Self::compute_tag(domain_separator, &io_pattern);

        let mut state = [F::ZERO; WIDTH];
        state[0] = tag;

        let rate = WIDTH - CAPACITY;

        let mut poseidon = Poseidon::new(curve);
        poseidon.state = state;

        Self {
            absorb_pos: 0,
            squeeze_pos: 0,
            io_count: 0,
            io_pattern,
            rate,
            capacity: CAPACITY,
            tag,
            poseidon,
        }
    }

    fn encode_io(io: &SpongeOp) -> u32 {
        match io {
            SpongeOp::Absorb(n) => (n + 0x80000000) as u32,
            SpongeOp::Squeeze(n) => (*n) as u32,
        }
    }

    // Compute tag as described in section 2.3 of the SAFE documentation
    pub fn compute_tag(domain_separator: &[u8], io_pattern: &IOPattern) -> F {
        // step 1: Encode
        let io_words = io_pattern
            .0
            .iter()
            .map(|io_i| Self::encode_io(io_i))
            .collect::<Vec<u32>>();

        // step 2: Aggregate
        let mut io_words_aggregated = vec![];
        for io_word in io_words {
            if io_words_aggregated.len() == 0 {
                io_words_aggregated.push(io_word);
            } else {
                let i = io_words_aggregated.len() - 1;
                if io_words_aggregated[i] > 0x80000000 && io_word > 0x80000000 {
                    io_words_aggregated[i] += io_word - 0x80000000;
                } else if io_words_aggregated[i] < 0x80000000 && io_word < 0x80000000 {
                    io_words_aggregated[i] += io_word;
                } else {
                    io_words_aggregated.push(io_word);
                }
            }
        }

        // step 3: Serialize
        let mut io_bytes = vec![];
        for io_word in io_words_aggregated {
            io_word.to_be_bytes().iter().for_each(|x| io_bytes.push(*x));
        }
        // Append the domain separator
        io_bytes.extend_from_slice(domain_separator);

        // step 4: Hash
        let mut hasher = Keccak::v256();
        hasher.update(&io_bytes);
        let mut result = [0u8; 32];
        hasher.finalize(&mut result);

        // The tag will be the first 128 bits of the hash.
        let mut tag = Vec::with_capacity(32);
        tag.extend_from_slice(&result[0..16]);

        // And we pad the tag to 256 bits to convert it to a field element.
        // TODO: Support variable field size
        tag.extend_from_slice(&[0; 16]);

        F::from_be_bytes_mod_order(&tag)
    }

    pub fn absorb(&mut self, x: &[F]) {
        if x.len() == 0 {
            return;
        }

        for x_i in x {
            if self.absorb_pos == self.rate {
                self.permute();
                self.absorb_pos = 0
            }

            self.poseidon.state[self.absorb_pos] = *x_i;
            self.absorb_pos += 1;
        }

        // assert_eq!(self.io_pattern.0[self.io_count], SpongeOp::Absorb(x.len()));

        self.io_count += 1;
        self.squeeze_pos = self.rate;
    }

    pub fn squeeze(&mut self, length: usize) -> Vec<F> {
        let mut y = Vec::with_capacity(length);
        if length == 0 {
            return vec![];
        }

        for _ in 0..length {
            if self.squeeze_pos == self.rate {
                self.permute();
                self.squeeze_pos = 0;
                self.absorb_pos = 0;
            }

            y.push(self.poseidon.state[self.squeeze_pos]);
            self.squeeze_pos += 1;
        }

        self.io_count += 1;
        y
    }

    pub fn finish(&self) -> Result<(), String> {
        if self.io_count != self.io_pattern.0.len() {
            return Err("IO pattern mismatch".to_string());
        }

        Ok(())
    }

    fn permute(&mut self) {
        self.poseidon.permute();
        self.poseidon.pos = 0;
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    type Fp = ark_secp256k1::Fq;
    const RATE: usize = 8;
    const WIDTH: usize = RATE + CAPACITY;

    #[test]
    fn test_poseidon_secp256k1() {
        let input = (0..RATE).map(|i| Fp::from(i as u64)).collect::<Vec<Fp>>();

        let mut sponge =
            PoseidonSponge::<Fp, WIDTH>::new(b"test", PoseidonCurve::SECP256K1, IOPattern(vec![]));
        sponge.absorb(&input);
        let digest = sponge.squeeze(1)[0];

        assert_eq!(
            digest.to_string(),
            "70319031943286297975812718474954003309712593879219035486413267041586554011143"
        );
    }

    #[test]
    fn test_poseidon_sponge() {
        let io_pattern = IOPattern(vec![
            SpongeOp::Absorb(2),
            SpongeOp::Squeeze(1),
            SpongeOp::Absorb(1),
            SpongeOp::Squeeze(3),
        ]);

        let inputs = vec![vec![Fp::from(1), Fp::from(2)], vec![Fp::from(3)]].concat();

        let mut sponge =
            PoseidonSponge::<Fp, WIDTH>::new(b"test", PoseidonCurve::SECP256K1, io_pattern.clone());

        let mut input_position = 0;
        for op in io_pattern.0 {
            match op {
                SpongeOp::Absorb(l) => {
                    sponge.absorb(&inputs[input_position..(input_position + l)]);
                    input_position += l;
                }
                SpongeOp::Squeeze(l) => {
                    sponge.squeeze(l);
                }
            }
        }

        assert_eq!(sponge.finish(), Ok(()));
    }
}
