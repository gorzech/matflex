% Compare calculated mass matrix with simple integration

% Define material mat and section properties sec
l = 0.5;
a = 0.05; % square cross section
sec = h_section(a);
mat = h_material();

n_elem = 16;
elem_len = l / n_elem;

[q0, elem_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);

%   fbody contains (at lest)
%   B2 - reference conditions
fbody.q0 = q0;
fbody.elem_idx = elem_idx;
fbody.mat = mat;
fbody.sec = sec;
fbody.elem_len = elem_len;
fbody.n_elem = n_elem;
fbody.mass = mat.rho * sec.A * l;
B2 = speye(size(M));
% Clamped floating frame at first node

% boundary conditions
% Simply supported beam - reasonably similar to free-free
bc_nodes = [1, 17];
bc_dofs = { 1:4, 2:3 };
[~, constr_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, length(q0));

B2(:, constr_dof) = [];
fbody.B2 = B2;

fbody.m_ff = B2' * M * B2;

%% Compare with simplified version - direct integration of L

% L = [I -A*us A*S]
p = rand(4, 1);
p = p ./ norm(p);
qf = rand(length(q0) - 6, 1) * 0.1;

M_ffr = ffr_mass_integrate(fbody, p, qf);
qn = q0 + B2 * qf;
M_expected = zeros(length(q0));
% Loop over elements!
for ii = 1 : n_elem
    mi = @(xi, eta, zeta) mass_int_L(xi, eta, zeta, B2, qn, p, elem_len, elem_idx(:, ii));
    M_expected = M_expected + volumetric_integral(mi, mat, sec, elem_len, 4, 2, 2);
end

diffM = M_expected - M_ffr;
M_ffr(4, 4) = 0; % This part may be needed Ixx
n_diffM2 = norm(M_expected - M_ffr);
assert(norm(diffM) < 1e-14 || n_diffM2 < 1e-14)

function M_int = mass_int_L(xi, eta, zeta, B2, qn, p, elem_len, elem_idx)
    S_elem = beam3d_shape_fun(xi, eta, zeta, elem_len);
    u = S_elem * qn(elem_idx);
    A = Rot(p);
    
    L = [eye(3), -A * sksym(u), A * S_elem * B2(elem_idx, :)];
    M_int = L' * L;
end
