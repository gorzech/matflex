% Test for quadratic inertia forces - simplified version. Even if version
% is reasonably simple, it is worth to verify an implementation

% Define material mat and section properties sec
mat.rho = 7801;
l = 0.5;
a = 0.05; % square cross section
sec.A = a * a;
sec.Iy = a ^ 4 / 12;
sec.Iz = sec.Iy;
sec.Jx = a ^ 4 / 6; % polar moment of inertia
sec.thky = a;
sec.thkz = a;

mat.E = 2e11;
nu = 0.3;
mat.G = mat.E / 2 / (1 + nu);

n_elem = 8;
elem_len = l / n_elem;
[q0, elem_idx, node_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);

B2 = speye(size(M));

% boundary conditions
% Simply supported beam
bc_nodes = [1, n_elem + 1];
bc_dofs = { 1:4, 2:3 };
[~, constr_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, length(q0));

B2(:, constr_dof) = [];

fbody = h_fbody2(q0, M, zeros(size(M)), B2, node_idx);
% fbody.m9 = [];

rng(9)
p = rand(4, 1); 
p = p ./ norm(p);
qf = rand(length(q0) - 6, 1) * 0.1;
% For zero velocity should be zero 
r_p = zeros(3, 1);
om = r_p;
qf_p = zeros(size(qf));
Qv0 = ffr_quadratic_inertia_preintegrated(fbody, p, qf, om, qf_p);

assert(all(size(Qv0) == [size(B2, 1), 1]), "Incorrect size of a Qv")
assert(norm(Qv0) < 1e-15, "Qv for zero velocity must be zero")

%% Compare with approximated version - same as first test in Qv simple integrate

% L = [I -A*us A*S]

v_coeff = 1.234;
r_p = rand(3, 1) .* v_coeff;
om = rand(3, 1) .* v_coeff;
qf_p = rand(size(qf)) .* 0.5;
q = [1; 1; 1; p; qf];

L = Gep(p);
p_p = 0.5 .* (L' * om);
q_p = [r_p; p_p; qf_p];

% Ensure transformation to work both ways
assert(norm(om - 2 * L * p_p) < 1e-15)

Qv1_2 = approxJacobian(@(iq) first_Qv(iq, fbody, r_p, p_p, qf_p), q, 1e-7);

Qv_expected1_A = - Qv1_2(1 : end - 1, :) * q_p;

Qv2 = 0.5 .* Qv1_2(end, :)';
Qv2_B = 0.5 .* Qv1_2(1 : end - 1, :)' * q_p;

% Two versions of getting second term of the Qv should result in the same
% vector
diff_Qv_2 = Qv2 - Qv2_B;
assert(norm(diff_Qv_2) < 1e-7, ...
    '%g - second Qv approx term from two versions of computations', norm(diff_Qv_2))

U = blkdiag(eye(3), 0.5 * L, eye(length(qf)));
% T = blkdiag(eye(3), 2 * L, eye(length(qf)));

% W.R.T q_p coordinates (that is Euler params)
Qv_expected_q = Qv_expected1_A + Qv2;

% Transformation to h coordinates - for angular velocity
Qv_expected = U * Qv_expected_q;

% Compute original vector
Qv = ffr_quadratic_inertia_preintegrated(fbody, p, qf, om, qf_p);

diffQv = Qv_expected - Qv;
err_Qv = norm(diffQv);
assert(err_Qv < 1e-6, "Too large error: %g", err_Qv)

function Qv1_2 = first_Qv(q, fbody, r_p, p_p, qf_p)
    p = q(4 : 7) ./ norm(q(4 : 7));
    L = Gep(p);
    T = blkdiag(eye(3), 2 * L, eye(length(qf_p)));
    M = ffr_mass_preintegrated(fbody, p, q(8 : end));
    M_t = T' * M * T;
    q_p = [r_p; p_p; qf_p];
    Qv_1 = M_t * q_p;
    Qv1_2 = [Qv_1
        q_p' * Qv_1];
end