use super::tree::CommittedMerkleTree;
use crate::FieldGC;

#[derive(Clone)]
pub struct TensorCode<F: FieldGC>(pub Vec<Vec<F>>);

impl<F: FieldGC> TensorCode<F> {
    pub fn commit(&self, num_cols: usize, num_rows: usize) -> CommittedTensorCode<F> {
        // Flatten the tensor codeword in column major order
        let mut tensor_codeword = vec![];
        for j in 0..(num_cols * 2) {
            for i in 0..num_rows {
                tensor_codeword.push(self.0[i][j])
            }
        }
        // Merkle commit the codewords
        let committed_tree = CommittedMerkleTree::from_leaves(tensor_codeword, num_cols * 2);

        CommittedTensorCode {
            committed_tree,
            tensor_codeword: Self(self.0.clone()),
        }
    }
}

#[derive(Clone)]
pub struct CommittedTensorCode<F: FieldGC> {
    pub committed_tree: CommittedMerkleTree<F>,
    pub tensor_codeword: TensorCode<F>,
}

impl<F: FieldGC> CommittedTensorCode<F> {
    pub fn query_column(&self, column: usize) -> Vec<F> {
        let num_rows = self.tensor_codeword.0.len();

        let leaves =
            self.committed_tree.leaves[column * num_rows..((column + 1) * num_rows)].to_vec();
        leaves
    }
}
