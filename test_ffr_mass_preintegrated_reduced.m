% Compare calculated mass matrix with simple integration for reduced system

% Define material mat and section properties sec
rng(3)
mat.rho = 7801;
l = 2;

a = 0.05; % square cross section
sec = h_section(a);

mat.E = 2e11;
nu = 0.3;
mat.G = mat.E / 2 / (1 + nu);

n_elem = 40;
elem_len = l / n_elem;
[q0, elem_idx, node_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);
K = fem_stiffness(elem_idx, mat, sec, elem_len);

n_modal = 2;
cb_nodes = [1, 41];
cb_dofs = {1:6, 1:6};
V = fem_CraigBampton_orth(n_modal, K, M, cb_nodes, cb_dofs);

% Apply mean axis conditions
V = V(:, 7:end);

for ii = 1 : 41
    q0(node_idx(1, ii)) = q0(node_idx(1, ii)) - 1;
end

fbody = h_fbody2(q0, M, K, V, node_idx);

%% Compare with simplified version - direct integration of L

% L = [I -A*us A*S]
p = rand(4, 1);
p = p ./ norm(p);
% p = [1; 0; 0; 0];
qf = rand(size(V, 2), 1) * 0.01;
% qf = rand(size(V, 2), 1) * 0.0;

M_ffr = full(ffr_mass_preintegrated(fbody, p, qf));
qn = q0 + V * qf;
M_expected = zeros(length(qf) + 6);
% Loop over elements!
for ii = 1 : n_elem
    mi = @(xi, eta, zeta) mass_int_L(xi, eta, zeta, V, qn, p, elem_len, elem_idx(:, ii));
    M_expected = M_expected + volumetric_integral(mi, mat, sec, elem_len, 4, 2, 2);
end
% Fix for the beam model
M_expected(4, 4) = M_expected(4, 4) + fbody.m7(1, 1);
% Assume that "exact" values is as for the beam - m * l ^ 2 / 12
M_expected(5, 5) = M_expected(5, 5) + 0.5 * fbody.m7(1, 1);
M_expected(6, 6) = M_expected(6, 6) + 0.5 * fbody.m7(1, 1);

% Translational part
diffM = M_expected - M_ffr;
diffM(4:end, 4:end) = 0;
assert(norm(diffM) < 2e-5)

% Rotational part - most differences
diffM = M_expected - M_ffr;
diffM(1:3, 1:3) = 0;
diffM(7:end, 7:end) = 0;
assert(norm(diffM) < 5e-3, "Error of the rotational part is %g", norm(diffM))

% Fully flexible part
diffM = M_expected(7:end, 7:end) - M_ffr(7:end, 7:end);

assert(norm(diffM) < 1e-15)


function M_int = mass_int_L(xi, eta, zeta, B2, qn, p, elem_len, elem_idx)
    S = zeros(3, size(B2, 1));
    S_elem = beam3d_shape_fun(xi, eta, zeta, elem_len);
    S(:, elem_idx) = S_elem;
    u = S * qn;
    A = Rot(p);
    
    L = [eye(3), -A * sksym(u), A * S * B2];
    M_int = L' * L;
end
