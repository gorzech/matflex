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

% boundary conditions
% Simply supported beam
bc_nodes = [1, n_elem + 1];
bc_dofs = { 1:4, 2:3 };
[~, constr_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, length(q0));

B2(:, constr_dof) = [];
fbody.B2 = B2;

fbody.m_ff = B2' * M * B2;

rng(2)
p = rand(4, 1);
p = p ./ norm(p);
qf = rand(length(q0) - 6, 1) * 0.1;
% For zero velocity should be zero 
r_p = zeros(3, 1);
om = r_p;
qf_p = zeros(size(qf));
Qv0 = ffr_quadratic_inertia_simple_integrate(fbody, p, qf, om, qf_p);

assert(all(size(Qv0) == [size(B2, 1), 1]), "Incorrect size of a Qv")
assert(norm(Qv0) < 1e-15, "Qv for zero velocity must be zero")

%% Compare with approximated version - using approx Jacobian

% Note that consistency is super important here. When working with
% different set of coordinates, they have to be employed up to the end of
% derivations!

% L = [I -A*us A*S]

v_coeff = 1.234;
r_p = rand(3, 1) .* v_coeff;
om = rand(3, 1) .* v_coeff;
qf_p = rand(size(qf)) .* v_coeff;
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
assert(norm(diff_Qv_2) < 2e-7, ...
    '%g - second Qv approx term from two versions of computations', norm(diff_Qv_2))

U = blkdiag(eye(3), 0.5 * L, eye(length(qf)));
T = blkdiag(eye(3), 2 * L, eye(length(qf)));

% W.R.T q_p coordinates (that is Euler params)
Qv_expected_q = Qv_expected1_A + Qv2;
% Needed only for comparison with T' *  Qv !
Qv_expected_q(4:7) = (eye(4) - p * p') * Qv_expected_q(4:7);

% Transformation to h coordinates - for angular velocity
% Note that this transformation nullify above transformation (as L * p = 0
Qv_expected = U * Qv_expected_q;

% Compute original vector
Qv = ffr_quadratic_inertia_simple_integrate(fbody, p, qf, om, qf_p);

% Compare vectors in q coordinates
diffQv_q = Qv_expected_q - T' *  Qv;
assert(norm(diffQv_q) < 2e-6, '%g - Difference in Qv in q (Euler params)', norm(diffQv_q))

diffQv = (Qv_expected - Qv);

err_Qv = norm(diffQv);
assert(err_Qv < 1e-6, "Too large error for Qv in h: %g", err_Qv)

%% Compare with integrated version - using time derivatives of L

v_coeff = 1.234;
r_p = rand(3, 1) .* v_coeff;
om = rand(3, 1) .* v_coeff;
qf_p = rand(size(qf)) .* v_coeff;
h = [r_p; om; qf_p];

LtL_p = zeros(length(q0));
for ii = 1 : n_elem
    qv = @(xi, eta, zeta) LtL_p_int_L(xi, eta, zeta, elem_idx(:, ii), ...
        B2, elem_len, q0 + B2 * qf, qf_p, Rot(p), sksym(om));
    LtL_p = LtL_p + volumetric_integral(qv, mat, sec, elem_len, 4, 2, 2);
end

Qv_expected1_B = - (LtL_p + LtL_p') * h; % Term -M_p * h
Qv_expected2_B = LtL_p' * h; % Term  + 0.5 * d/d_q (h * M * h);

Qv = ffr_quadratic_inertia_simple_integrate(fbody, p, qf, om, qf_p);

diffQv = Qv_expected1_B + Qv_expected2_B - Qv;

err_Qv = norm(diffQv);
assert(err_Qv < 1e-14, "Too large error: %g", err_Qv)

function Qv1_2 = first_Qv(q, fbody, r_p, p_p, qf_p)
    p = q(4 : 7) ./ norm(q(4 : 7));
    L = Gep(p);
    T = blkdiag(eye(3), 2 * L, eye(length(qf_p)));
    M = ffr_mass_integrate(fbody, p, q(8 : end));
    M_t = T' * M * T;
    q_p = [r_p; p_p; qf_p];
    Qv_1 = M_t * q_p;
    Qv1_2 = [Qv_1
        q_p' * Qv_1];
end

function LtL_p_int = LtL_p_int_L(xi, eta, zeta, elem_idx, B2, elem_len, qn, qf_p, A, oms)
    S_elem = beam3d_shape_fun(xi, eta, zeta, elem_len);
    
    u = S_elem * qn(elem_idx);
    us = sksym(u);
    
    phi = S_elem * B2(elem_idx, :);
    u_p = phi * qf_p;
    us_p = sksym(u_p);
    
    L = [eye(3), -A * us, A * phi];
    L_p = [zeros(3), -A * (oms * us + us_p), A * oms * phi];
    LtL_p_int = L' * L_p;
end