function [qf0, K, M, node_idx] = h_fem_model(body_len, n_elem, sec, mat, ...
    lumped_mass)
%H_FEM_MODEL Helper to get beam FEM model matrices

p = inputParser;
p.StructExpand = false;
positive_scalar = @(x) isscalar(x) && isnumeric(x) && x > 0;
addRequired(p, 'body_len', positive_scalar);
addOptional(p, 'n_elem', 2, positive_scalar);
addOptional(p, 'sec', h_section());
addOptional(p, 'mat', h_material());
addOptional(p, 'lumped_mass', false);
parse(p, body_len, n_elem, sec, mat, ...
    lumped_mass);

elem_len = body_len / p.Results.n_elem;
[qf0, elem_idx, node_idx] = straight_beam3d(elem_len, p.Results.n_elem);

M = fem_mass(elem_idx, p.Results.mat, p.Results.sec, elem_len, p.Results.lumped_mass);
K = fem_stiffness(elem_idx, p.Results.mat, p.Results.sec, elem_len);

end

