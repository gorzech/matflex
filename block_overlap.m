function B = block_overlap(block_idx, B_block)
%BLOCK_OVERLAP Generate block overlap sprase matrix like required e.g., in
%FEM analysis
%   block_idx - column has indices of the block in B
%   B_block - dense matrix to be blocked
%   B - output sparse matrix

rows = size(B_block, 1);
blocks = size(block_idx, 2);

% generate indices for sparse matrix as required for sparse matrix
i_idx = repmat(block_idx, [rows, 1]);
j_idx = repmat(block_idx(:)', [rows, 1]);
B_flat = repmat(B_block, [1, blocks]);

B = sparse(i_idx(:), j_idx(:), B_flat(:));