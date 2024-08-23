function [V, Kv, freq] = fem_CraigBampton_orth(n_modal, K, M, cb_nodes, cb_dofs)
%FEM_CRAIGBAMPTON_ORTH Computes transformation matrix for CraigBampton
%method. This is done WITH orthonormalization.
%   n_modal - number of modal modes to be added to static modes
%   K, M - system stiffness and mass matrices
%   cb_nodes, cb_dofs - define CB boundary dofs as required by fem_bc_split
%   Output
%   V - Craig-Bampton transformation matrix
%   Kv - optional reduced system matrix Mv should be unity
%   freq - optional natural frequencies

% First stage CB
[V_I, Kv_I, Mv_I] = fem_CraigBampton(n_modal, K, M, cb_nodes, cb_dofs);

% orthonormalize - makes them symmetric as they may loose this property
[eigV, freq] = normalized_eig(0.5 * (Kv_I' + Kv_I), 0.5 * (Mv_I' + Mv_I)); 

V = V_I * eigV;

if nargout > 1
    Kv = eigV' * Kv_I * eigV;
end

end

