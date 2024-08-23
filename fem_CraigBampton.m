function [V, Kv, Mv, freq] = fem_CraigBampton(n_modal, K, M, cb_nodes, cb_dofs)
%FEM_CRAIGBAMPTON Computes transformation matrix for CraigBampton
%method. This is done without orthonormalization.
%   n_modal - number of modal modes to be added to static modes
%   K, M - system stiffness and mass matrices
%   cb_nodes, cb_dofs - define CB boundary dofs as required by fem_bc_split
%   Output
%   V - Craig-Bampton transformation matrix
%   Kv, Mv - optional reduced system matrices

[K_ii, ~, K_ib, ~, P] = fem_bc_split(cb_nodes, cb_dofs, K);
M_ii = fem_bc_split(cb_nodes, cb_dofs, M);

% Static correction modes - apply unit force at boundary
Vc = -K_ii \ K_ib;

% For 'chol' option eigenvectors are (usually) M-normalized
[eigV, freq] = normalized_eig(K_ii, M_ii, n_modal); 

n_dofs = size(K, 1);
n_bnd = size(K_ib, 2);

V = zeros(n_dofs, n_bnd + n_modal);
V(P, 1:n_bnd) = [eye(n_bnd); Vc];
V(P(n_bnd + 1 : end), n_bnd + 1 : end) = eigV;

if nargout > 1
    Kv = V' * K * V;
    if nargout > 2
        Mv = V' * M * V;
    end 
end

end

