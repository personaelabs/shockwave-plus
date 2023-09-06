use super::utils::hash_all;
use crate::FieldGC;
use ark_ff::BigInteger;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[derive(Clone)]
pub struct CommittedMerkleTree<F> {
    pub column_roots: Vec<[u8; 32]>,
    pub leaves: Vec<F>,
    pub num_cols: usize,
    pub root: [u8; 32],
}

impl<F: FieldGC> CommittedMerkleTree<F> {
    pub fn from_leaves(leaves: Vec<F>, num_cols: usize) -> Self {
        let n = leaves.len();
        let num_rows = n / num_cols;

        #[cfg(not(feature = "parallel"))]
        let leaf_bytes = leaves
            .iter()
            .map(|x| x.into_bigint().to_bytes_be().try_into().unwrap())
            .collect::<Vec<[u8; 32]>>();

        #[cfg(not(feature = "parallel"))]
        let column_roots = (0..num_cols)
            .map(|col| {
                let column_leaves = leaf_bytes[col * num_rows..(col + 1) * num_rows].to_vec();
                let column_root = hash_all(&column_leaves);
                column_root
            })
            .collect::<Vec<[u8; 32]>>();

        #[cfg(feature = "parallel")]
        let leaf_bytes = leaves
            .par_iter()
            .map(|x| x.into_bigint().to_bytes_be().try_into().unwrap())
            .collect::<Vec<[u8; 32]>>();

        #[cfg(feature = "parallel")]
        let column_roots = (0..num_cols)
            .into_par_iter()
            .map(|col| {
                let column_leaves = leaf_bytes[col * num_rows..(col + 1) * num_rows].to_vec();
                let column_root = hash_all(&column_leaves);
                column_root
            })
            .collect::<Vec<[u8; 32]>>();

        let root = hash_all(&column_roots);

        Self {
            column_roots,
            leaves,
            root,
            num_cols,
        }
    }

    pub fn root(&self) -> [u8; 32] {
        self.root
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct BaseOpening {
    pub hashes: Vec<[u8; 32]>,
}

impl BaseOpening {
    pub fn verify(&self, root: [u8; 32]) -> bool {
        let r = hash_all(&self.hashes);
        root == r
    }
}
