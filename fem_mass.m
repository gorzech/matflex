function M = fem_mass(elem_idx, mat, sec, elem_length, lumped)
%FEM_MASS Generate system mass for a straight beam
%   Input:
%   elem_idx - indices of the elements in system coordinate vector (from
%   e.g. straight_beam script)
%   mat - material properties
%   sec - section properties (area, moments of inertia, etc.)
%   elem_length - length of a single element - uniform mesh
%   lumped - use lumped mass approx (default to false)
%   Output:
%   M - system mass matrix. Sparse

if nargin < 5
    lumped = false;
end

% As we assume uniform beam mesh, only one mass matrix evaluation is
% required
M_elem = beam3d_mass(mat.rho, sec.A, elem_length, sec.Jx, sec.Iy, sec.Iz, lumped);

M = block_overlap(elem_idx, M_elem);

end

