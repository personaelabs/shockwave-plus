use core::num;
use std::marker::PhantomData;

use super::utils::hash_two;
use crate::{utils::hash_all, FieldExt};
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct CommittedMerkleTree<F> {
    pub column_roots: Vec<[u8; 32]>,
    pub leaves: Vec<F>,
    pub num_cols: usize,
    pub root: [u8; 32],
}

impl<F: FieldExt> CommittedMerkleTree<F> {
    pub fn from_leaves(leaves: Vec<F>, num_cols: usize) -> Self {
        let n = leaves.len();
        debug_assert!(n.is_power_of_two());
        let num_rows = n / num_cols;
        assert!(num_rows & 1 == 0); // Number of rows must be even

        let leaf_bytes = leaves
            .iter()
            .map(|x| x.to_repr())
            .collect::<Vec<[u8; 32]>>();

        let mut column_roots = Vec::with_capacity(num_cols);
        for col in 0..num_cols {
            let column_leaves = leaf_bytes[col * num_rows..(col + 1) * num_rows].to_vec();
            let column_root = hash_all(&column_leaves);
            column_roots.push(column_root);
        }

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

    pub fn leaves(&self) -> Vec<F> {
        self.leaves.clone()
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct BaseOpening {
    pub hashes: Vec<[u8; 32]>,
}

impl BaseOpening {
    pub fn verify(&self, root: [u8; 32]) -> bool {
        let r = hash_all(&self.hashes);
        root == r
    }
}
