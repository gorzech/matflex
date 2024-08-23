function [K, M, Fg, qf, node_idx] = th_beam_fem_model(grav, lumped)
% 2 m flexible beam - for different test purposes

if nargin < 1
    grav = zero(3, 1);
end
if nargin < 2
    lumped = false;
end
% Define material mat and section properties sec
[sec, mat, l, n_elem] = th_beam_fem_props();

elem_len = l / n_elem;
[~, elem_idx] = straight_beam3d(elem_len, n_elem);

[qf, K, M, node_idx] = h_fem_model(l, n_elem, ...
    sec, mat, lumped);

if nargout > 2
    Fg_element = beam3d_grav(mat.rho, sec.A, elem_len, grav, lumped);
    Fg = block_vector(elem_idx, Fg_element);
end