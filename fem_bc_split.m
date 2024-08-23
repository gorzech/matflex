function [M_ii, M_bb, M_ib, M_bi, P] = fem_bc_split(bc_nodes, bc_dofs, M)
%FEM_BC_SPLIT Apply boundary conditions to a system matrices.
%Works for my beam only. Splot according to submatrices:
%               /            \
%               | M_bb  M_bi |
%    M(P, P)  = |            |
%               | M_ib  M_ii |
%               \            /
%       and b indicate defined boundary conditions and P is permutation
%       vector
%   Input:
%   bc_nodes - list of nodes where to apply conditions
%   bc_dofs - cell array of size of node_id with associated dofs to be removed
%   (1 to 6 for our beam)
%   M - matrix to be transformed
%   Output:
%   M_ii - submatrix with removed DOFs
%   M_bb - submatrix with preserved boundary DOFs
%   M_ib - submatrix in lower left part
%   M_bi - submatrix in upper right part
%   P - permutation vector - boundary dofs first, next internal dofs

[internal_dofs, boundary_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, size(M, 1));
% Simplest approach
M_ii = M(internal_dofs, internal_dofs);
if nargout > 1
    M_bb = M(boundary_dof, boundary_dof);
end
if nargout > 2
    M_ib = M(internal_dofs, boundary_dof);
end
if nargout > 3
    M_bi = M(boundary_dof, internal_dofs);
end
if nargout > 4
    P = [boundary_dof; internal_dofs];
end

end

