% Basic test for translational spring damper

%% Two rigid bodies

rng(104);

k = 1e4; c = 1e2;

r1 = [0; 1; 0];
p1 = normunit(rand(4, 1));
r2 = [0.1; 0.2; 0.7];
p2 = normunit(rand(4, 1));

sys = sys_generate([h_body(1, eye(3)), h_body(1, eye(3))], [], ...
    'usequadratic', false);

q0 = [r1; p1; r2; p2];

point1 = [1; 0; 0];
point2 = [2; 0; 0];

sys = h_sys_force_tsd(sys, k, c, 1, point1, point2, q0, 2, 1);

assert(length(sys.force.tsd) == 1);

Q = mbd_force(sys, 0.1, q0, zeros(sys.nh, 1));
assert(length(Q) == 12);
assert(norm(Q) < 2e-11)

q = q0;
q(1) = q(1) - 0.1;
h = zeros(sys.nh, 1);
h(1) = 1;

Q = mbd_force(sys, 0.1, q, h);

F = (k * 0.1 - c * 1) * [1; 0; 0];

Q_expected = [F;
    Rot(p1)' * sksym(point1 - [0.1; 0; 0] - r1) * F
    -F
    -Rot(p2)' * sksym(point2 - r2) * F];

assert(norm(Q - Q_expected) < 1e-10)

% Add extra force
sys = h_sys_force_tsd(sys, k, 0, 1, point1, point2, q, 1, 2);

assert(length(sys.force.tsd) == 2);
Q = mbd_force(sys, 0.1, q, h);
assert(norm(Q - Q_expected) < 1e-10)

%% Two flexible bodies

k = 17e4; d = 123.7;

rng(71);
p1 = normunit(rand(4, 1));
p2 = normunit(rand(4, 1));

fb1 = h_fbody(6, 4);
fb2 = h_fbody(0.2, 3);
sys = sys_generate([], [fb1, fb2]);

s1 = [9; 2; -1];

s2 = [10; 2; -1];

s1_local = [4.5; 0; 0];
s2_local = [fb2.q0(7); 0; 0];
r1 = s1 - Rot(p1) * s1_local;
r2 = s2 - Rot(p2) * s2_local;

q0 = [r1; p1; zeros(fb1.nflex, 1); r2; p2; zeros(fb2.nflex, 1)];

sys = h_sys_force_tsd(sys, k, d, 1, s1, s2, q0, 1, 2);

assert(length(sys.force.tsd) == 1);

Q = mbd_force(sys, 0.1, q0, zeros(sys.nh, 1));
assert(length(Q) == sys.nh);
assert(norm(Q) < 2e-14)