% Compare calculated mass matrix with simple integration

% Define material mat and section properties sec
rng(3)
l = 0.5;
a = 0.05; % square cross section
sec = h_section(a);
mat = h_material();

n_elem = 16;
elem_len = l / n_elem;

[q0, elem_idx, node_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);

B2 = speye(size(M));
% Clamped floating frame at first node

% boundary conditions
% Simply supported beam - reasonably similar to free-free
bc_nodes = [1, 17];
bc_dofs = { 1:4, 2:3 };
[~, constr_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, length(q0));
B2(:, constr_dof) = [];

fbody = h_fbody2(q0, M, zeros(size(M)), B2, node_idx);

%% Compare with simplified version - direct integration of L

% L = [I -A*us A*S]
p = rand(4, 1);
p = p ./ norm(p);
qf = rand(length(q0) - 6, 1) * 0.01;

M_ffr = ffr_mass_preintegrated(fbody, p, qf);
qn = q0 + B2 * qf;
M_expected = zeros(length(q0));
% Loop over elements!
for ii = 1 : n_elem
    mi = @(xi, eta, zeta) mass_int_L(xi, eta, zeta, B2, qn, p, elem_len, elem_idx(:, ii));
    M_expected = M_expected + volumetric_integral(mi, mat, sec, elem_len, 4, 2, 2);
end
% Fix for the beam model
M_expected(4, 4) = M_expected(4, 4) + fbody.m7(1, 1);

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
