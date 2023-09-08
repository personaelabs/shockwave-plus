use crate::tensor_pcs::hasher::Hasher;
use crate::FieldGC;

#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

#[derive(Clone)]
pub struct CommittedMerkleTree<F: FieldGC, H: Hasher<F>> {
    pub column_roots: Vec<H::T>,
    pub leaves: Vec<F>,
    pub num_cols: usize,
    pub root: H::T,
}

impl<F: FieldGC, H: Hasher<F>> CommittedMerkleTree<F, H> {
    pub fn from_leaves(leaves: Vec<F>, num_cols: usize, hasher: &H) -> Self {
        let n = leaves.len();
        let num_rows = n / num_cols;

        #[cfg(not(feature = "parallel"))]
        let column_roots = (0..num_cols)
            .map(|col| {
                let column_leaves = leaves[col * num_rows..(col + 1) * num_rows].to_vec();
                let column_root = hasher.hash_felts(&column_leaves);
                column_root
            })
            .collect::<Vec<[u8; 32]>>();

        #[cfg(feature = "parallel")]
        let column_roots = (0..num_cols)
            .into_par_iter()
            .map(|col| {
                let column_leaves = leaves[col * num_rows..(col + 1) * num_rows].to_vec();
                let column_root = hasher.hash_felts(&column_leaves);
                column_root
            })
            .collect::<Vec<H::T>>();

        let root = hasher.hash_all(&column_roots);

        Self {
            column_roots,
            leaves,
            root,
            num_cols,
        }
    }

    pub fn root(&self) -> H::T {
        self.root
    }
}
