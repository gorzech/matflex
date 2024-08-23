% Script used to test constraints for perp1 constraint

%% Start with two bodies

rng(8);
r1 = [0; 1; 0];
p1 = normunit(rand(4, 1));
r2 = [0.1; 0.2; 0.7];
p2 = normunit(rand(4, 1));

sys = sys_generate([h_body(1, eye(3)), h_body(1, eye(3))], []);

q0 = [r1; p1; r2; p2];

v_1 = [1; 0; 0];
v_2 = [0; 1; 0];

sys = h_sys_constr_perp1(sys, q0, 1, 2, v_1, v_2, []);

assert(sys.nconstr == 1);
assert(length(sys.constr.perp1) == 1);

[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 1);

assert(norm(c) < 1e-14);

%Translate and rotate (around 0) both bodies

rt = [17; -11.2; pi/4];
% rt = [0; 0; 0];
pt = q_axis(1.1324, 3);
% pt = q_axis(0, 3);
At = Rot(pt);

q = [At * (r1 + rt); q_mult(pt, p1); At * (r2 + rt); q_mult(pt, p2)];

[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c) < 1e-14);

% Now rotate and translate first body only
q = [At * (r1 + rt); q_mult(pt, p1); r2; p2];

% New vector after rotation:
v_1_n = At * v_1;

c_expected = dot(v_1_n, v_2);
[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c - c_expected) < 1e-14);

%% Two flexible bodies

rng(8);
p1 = normunit(rand(4, 1));
p2 = normunit(rand(4, 1));

fb1 = h_fbody(6, 4);
fb2 = h_fbody(0.2, 3);
sys = sys_generate([], [fb1, fb2]);

s0 = [9; 2; -1];

s1_local = [4.5; 0; 0];
s2_local = [fb2.q0(7); 0; 0];
r1 = s0 - Rot(p1) * s1_local;
r2 = s0 - Rot(p2) * s2_local;

q0 = [r1; p1; zeros(fb1.nflex, 1); r2; p2; zeros(fb2.nflex, 1)];

v_1 = [1; 0; 0];
v_2 = [0; 1; 0];
% For flexible bodies locations are needed
sys = h_sys_constr_perp1(sys, q0, 1, 2, v_1, v_2, s0);

assert(sys.nconstr == 1);
assert(length(sys.constr.perp1) == 1);

[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 1);

assert(norm(c) < 1e-14);

%Translate and rotate (around 0) both bodies

rt = [17; -11.2; pi/4];
pt = q_axis(1.1324, 3);
At = Rot(pt);

q = [At * (r1 + rt); q_mult(pt, p1)
    zeros(fb1.nflex, 1)
    At * (r2 + rt); q_mult(pt, p2)
    zeros(fb2.nflex, 1)];

[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c) < 1e-14);

% Now rotate and translate second body
q = [r1; p1
    zeros(fb1.nflex, 1)
    At * (r2 + rt); q_mult(pt, p2)
    zeros(fb2.nflex, 1)];

v_2_n = At * v_2;

c_expected = dot(v_1, v_2_n);
[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c - c_expected) < 1e-14);

% Add some deformation to the flexible body (first)
q = q0;
q(23:25) = 0.01;

At = eye(3) + sksym(q(23:25));
% Unit vector at the same body
A1 = Rot(p1);
v_1_n = A1 * At * A1' * v_1;

c_expected = dot(v_1_n, v_2);
[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c - c_expected) < 1e-14);

%% 4 bodies, 3 constraints - verify against approximation

rng(8);

fb1 = h_fbody(6, 4);
fb2 = h_fbody(0.2, 3);
sys = sys_generate([h_body(1, eye(3)), h_body(1, eye(3))], [fb1, fb2]);

r1_r = [0; 1; 0];
p1_r = normunit(rand(4, 1));
r2_r = [0.1; 0.2; 0.7];
p2_r = normunit(rand(4, 1));

s1_local = [4.5; 0; 0];
s2_local = [fb2.q0(7); 0; 0];
p1_f = normunit(rand(4, 1));
p2_f = normunit(rand(4, 1));

s0 = [-9; 12; 1/3];
r1_f = s0 - Rot(p1_f) * s1_local;
r2_f = s0 - Rot(p2_f) * s2_local;


q0 = [r1_r; p1_r; r2_r; p2_r
    r1_f; p1_f; zeros(fb1.nflex, 1); r2_f; p2_f; zeros(fb2.nflex, 1)];

v_1 = [1; 0; 0];
v_2 = [0; 0.3; 0.4];
% For flexible bodies locations are needed
sys = h_sys_constr_perp1(sys, q0, 1, 2, v_1, v_2, []);
sys = h_sys_constr_perp1(sys, q0, 4, 2, v_1, [0; -0.1; -0.07], s0);
sys = h_sys_constr_perp1(sys, q0, 3, 4, [0; 7; 3], [0; -3; 7], s0);

% Constr should be zero now
[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 3);
assert(norm(c) < 1e-14);

h = rand(sys.nh, 1) * 0.7623;
q0 = q0 + rand(sys.nq, 1) * 1e-2;
% normalize quaternions
for b = sys.bodies
    q0(b.q_idx(4 : 7)) = normunit(b.q_idx(4 : 7));
end
for b = sys.fbodies
    q0(b.q_idx(4 : 7)) = normunit(b.q_idx(4 : 7));
end

[~, c_q, c_p, g] = mbd_constr(sys, 0, q0, h);
[expected_c_h, expected_c_p, expected_g] = th_approx_constr_derivatives(sys, @mbd_constr, q0, h);

assert(norm(c_q - expected_c_h) / norm(c_q) < 4e-6)
assert(norm(c_p - expected_c_p) < 1e-14)
assert(norm(g - expected_g) < 4e-6)