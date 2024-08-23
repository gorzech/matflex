% Script to make basic test for drive point contraint

%% One rigid body

rng(29);
r1 = [0.1; 0.2; 0.7];
p1 = normunit(rand(4, 1));

sys = sys_generate(h_body(1, eye(3)), []);

q0 = [r1; p1];

s0 = r1 + [1.2; 0; -0.2];
c_idx = 1:2;
sys = h_sys_constr_drive_point(sys, s0, q0, 1, c_idx, ...
    @(t) [1.3; 0.2], @(t) [0; 0], @(t) [0; 0]);

assert(sys.nconstr == 2);
assert(length(sys.constr.drive_point) == 1);

[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 2);

assert(norm(c) < 1e-14);

% Translate and rotate body

rt = [17; -11.2; pi/4];
% rt = [0; 0; 0];
pt = q_axis(1.1324, 3);
% pt = q_axis(0, 3);
At = Rot(pt);

q = [At * (r1 + rt); q_mult(pt, p1)];

% Location of the joint point for the body:
st = At * (s0 + rt);

c_expected = st - s0;
[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c - c_expected(c_idx)) < 1e-14);

%% Flexible body

rng(665);
p1 = normunit(rand(4, 1));

fb1 = h_fbody(6, 4);
sys = sys_generate([], fb1);

s0 = [9; 2; -1];

s1_local = [4.5; 0; 0];
r1 = s0 - Rot(p1) * s1_local;

q0 = [r1; p1; zeros(fb1.nflex, 1)];

c_idx = [1, 3];
sys = h_sys_constr_drive_point(sys, s0, q0, 1, c_idx, @(t) s0(c_idx), ...
    @(t)[0; 0], @(t)[0; 0]);

assert(sys.nconstr == 2);
assert(length(sys.constr.drive_point) == 1);

[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 2);

assert(norm(c) < 1e-14);

% Translate and rotate (around 0) body

rt = [17; -11.2; pi/4];
% rt = [0; 0; 0];
pt = q_axis(1.1324, 3);
% pt = q_axis(0, 3);
At = Rot(pt);

q = [At * (r1 + rt); q_mult(pt, p1)
    zeros(fb1.nflex, 1)];

% Location of the joint point for the second body:
st = At * (s0 + rt);

c_expected = st - s0;
[c, ~, ~, ~] = mbd_constr(sys, 0, q, zeros(sys.nh, 1));
assert(norm(c - c_expected(c_idx)) < 1e-14);

%% 4 bodies, 4 constraints - verify against approximation
rng(965);

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

s0_3 = [-9; 12; 1/3];
r1_f = s0_3 - Rot(p1_f) * s1_local;
s0_4 = rand(3, 1);
r2_f = s0_4 - Rot(p2_f) * s2_local;


q0 = [r1_r; p1_r; r2_r; p2_r
    r1_f; p1_f; zeros(fb1.nflex, 1); r2_f; p2_f; zeros(fb2.nflex, 1)];

sys = h_sys_constr_drive_point(sys, s0_3, q0, 3, 1:2, ...
    @(t)s0_3(1:2) - t^3, @(t) [-3; -3] * t * t, ...
    @(t) [ -6; -6] * t);
sys = h_sys_constr_drive_point(sys, s0_4, q0, 4, 3, ...
    @(t)s0_4(3) + sin(t), @(t) cos(t), @(t) -sin(t));
s0_1 = rand(3, 1);
sys = h_sys_constr_drive_point(sys, s0_1, q0, 1, 1:3, ...
    @(t) s0_1 + (exp(t) - 1), @(t) ones(3, 1) .* exp(t), @(t) ones(3, 1) .* exp(t));
s0_2 = rand(3, 1);
sys = h_sys_constr_drive_point(sys, s0_2, q0, 1, [1, 3], ...
    @(t) s0_2([1, 3]), @(t) zeros(2, 1), @(t) zeros(2, 1));

% Constr should be zero now
[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 8);
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

tbase = 0.1;
[~, c_q, c_p, g] = mbd_constr(sys, tbase, q0, h);
[expected_c_h, expected_c_p, expected_g] = ...
    th_approx_constr_derivatives(sys, @mbd_constr, q0, h, tbase);

assert(norm(c_q - expected_c_h) / norm(c_q) < 4e-6)
assert(norm(c_p - expected_c_p) < 3.5e-7) % due to time approximation
assert(norm(g - expected_g) ./ norm(g) < 1e-6)