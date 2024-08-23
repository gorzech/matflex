% function qv
% Test for quadratic inertia forces

file_name = 'test_import_straight_beam.SID_FEM';
fbody = h_fbody_sid_fem(file_name);

rng(9)
p = rand(4, 1); 
p = p ./ norm(p);
qf = rand(fbody.nflex, 1) * 0.01;
% For zero velocity should be zero 
r_p = zeros(3, 1);
om = r_p;
qf_p = zeros(size(qf));
Qv0 = ffr_quadratic_inertia_sid(fbody, p, qf, om, qf_p);

assert(isequal(size(Qv0), [fbody.nh, 1]), "Incorrect size of a Qv")
assert(norm(Qv0) < 1e-15, "Qv for zero velocity must be zero")

%% Compare with approximated version - same as first test in Qv simple integrate

% L = [I -A*us A*S]

v_coeff = 1.234;
r_p = rand(3, 1) .* v_coeff;
om = rand(3, 1) .* v_coeff;
qf_p = rand(size(qf)) .* 0.1;
q = [1; 1; 1; p; qf];

L = Gep(p);
p_p = 0.5 .* (L' * om);
q_p = [r_p; p_p; qf_p];

% Ensure transformation to work both ways
assert(norm(om - 2 * L * p_p) < 1e-15)

Qv1_2 = approxJacobian(@(iq) first_Qv(iq, fbody, r_p, p_p, qf_p), q, 1e-6);

Qv_expected1_A = - Qv1_2(1 : end - 1, :) * q_p;

Qv2 = 0.5 .* Qv1_2(end, :)';
Qv2_B = 0.5 .* Qv1_2(1 : end - 1, :)' * q_p;

% Two versions of getting second term of the Qv should result in the same
% vector
diff_Qv_2 = Qv2 - Qv2_B;
assert(norm(diff_Qv_2) < 1e-7, ...
    '%g - second Qv approx term from two versions of computations', norm(diff_Qv_2))

U = blkdiag(eye(3), 0.5 * L, eye(length(qf)));

% W.R.T q_p coordinates (that is Euler params)
Qv_expected_q = Qv_expected1_A + Qv2;

% Transformation to h coordinates - for angular velocity
Qv_expected = U * Qv_expected_q;

% Compute original vector
Qv = ffr_quadratic_inertia_sid(fbody, p, qf, om, qf_p);

diffQv = Qv_expected - Qv;
% norm(diffQv(4 : 6))
% norm(diffQv(7 : end))
err_Qv = norm(diffQv);
assert(err_Qv < 3e-2, "Too large error: %g", err_Qv)
% end
function Qv1_2 = first_Qv(q, fbody, r_p, p_p, qf_p)
    p = q(4 : 7) ./ norm(q(4 : 7));
    L = Gep(p);
    T = blkdiag(eye(3), 2 * L, eye(length(qf_p)));
%     fbody.Cr = fbody.Cr.m0;
%     fbody.mmi = fbody.mmi.m0;
    M = ffr_mass_sid(fbody, p, q(8 : end));
    M_t = T' * M * T;
    q_p = [r_p; p_p; qf_p];
    Qv_1 = M_t * q_p;
    Qv1_2 = [Qv_1
        q_p' * Qv_1];
end