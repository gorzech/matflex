function Bv = block_vector(block_idx, v)
%BLOCK_OVERLAP Generate block overlap vector like required e.g., in
%FEM analysis
%   block_idx - column has indices of the block in Bv
%   v - vector to be blocked
%   Bv - output vector of larger size

n_elem = size(block_idx, 2);
max_idx = max(max(block_idx));

Bv = zeros(max_idx, 1);

for ii = 1 : n_elem
    Bv(block_idx(:, ii)) = Bv(block_idx(:, ii)) + v;
end