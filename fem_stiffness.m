function K = fem_stiffness(elem_idx, mat, sec, elem_length)
%FEM_STIFFNESS Stiffness matrix for a system with straight beam
%   Input:
%   elem_idx - indices of the elements in system coordinate vector (from
%   e.g. straight_beam script)
%   mat - material properties
%   sec - section properties (area, moments of inertia, etc.)
%   elem_length - length of a single element - uniform mesh
%   Output:
%   K - system stiffness matrix. Sparse

K_elem = beam3d_stiffness(mat.E, mat.G, sec.A, elem_length, sec.Jx, sec.Iy, sec.Iz);

K = block_overlap(elem_idx, K_elem);
end

