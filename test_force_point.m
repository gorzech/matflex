% Have to add some bodies to get forces

%% Rigid body

rng(12);
r1 = [0; 1; 0];
p1 = normunit(rand(4, 1));

sys = sys_generate(h_body(1, eye(3)), [], 'usequadratic', false);

q0 = [r1; p1];

force_point = [0; 1; -1];
s0 = r1 + [0; 1; -1];

force = @(t) [3; 2 + 2*t; t^2 - 1];

sys = h_sys_force_point(sys, force, s0, q0, 1);

assert(length(sys.force.point) == 1);

Q = mbd_force(sys, 0.1, q0, zeros(sys.nh, 1));
assert(length(Q) == 6);

F = force(0.1);
Q_expected = [F;
    Rot(p1)' * sksym(force_point) * F];

assert(norm(Q - Q_expected) < 1e-14)

% Add extra force

sys = h_sys_force_point(sys, @(t)-force(t), r1, q0, 1);

assert(length(sys.force.point) == 2);

Q_expected(1 : 3) = 0;
Q = mbd_force(sys, 0.1, q0, zeros(sys.nh, 1));
assert(norm(Q - Q_expected) < 1e-14)

%% 2 bodies

rng(8);

fb1 = h_fbody(6, 4);
sys = sys_generate(h_body(1, eye(3)), fb1, 'usequadratic', false);

r1_r = [0; 1; 0];
p1_r = normunit(rand(4, 1));

s1_local = [4.5; 0; 0]; 
p1_f = normunit(rand(4, 1));

s0 = [-9; 12; 1/3];
r1_f = s0 - Rot(p1_f) * s1_local;

q0 = [r1_r; p1_r
    r1_f; p1_f; zeros(fb1.nflex, 1)];

force = @(t) [-11; 125*t + pi; 2 - t ^ 3];

sys = h_sys_force_point(sys, force, s0, q0, 2);

F = force(0.17);
Q = mbd_force(sys, 0.17, q0, zeros(sys.nh, 1));

Q_expected = [zeros(6, 1); 
    F
    Rot(p1_f)' * sksym(s0 - r1_f) * F
    fb1.B2(19:21, :)' * Rot(p1_f)' * F];

assert(norm(Q - Q_expected) < 1e-13)

%% SID flexible body + rigid body

rng(8);

file_name = 'test_import_straight_beam.SID_FEM';
fb1 = h_fbody_sid_fem(file_name);
sys = sys_generate(h_body(1, eye(3)), fb1, 'usequadratic', false);

r1_r = [0; 1; 0];
p1_r = normunit(rand(4, 1));


s1_local = [-1; 0; 0]; 
sid_node_id = 1; % For that local point
p1_f = normunit(rand(4, 1));

s0 = [-9; 12; 1/3];
r1_f = s0 - Rot(p1_f) * s1_local;

q0 = [r1_r; p1_r
    r1_f; p1_f; zeros(fb1.nflex, 1)];

force = @(t) [-11; 125*t + pi; 2 - t ^ 3];

sys = h_sys_force_point(sys, force, s0, q0, 2);

F = force(0.17);
Q = mbd_force(sys, 0.17, q0, zeros(sys.nh, 1));

Q_expected = [zeros(6, 1); 
    F
    Rot(p1_f)' * sksym(s0 - r1_f) * F
    fb1.node(sid_node_id).origin.m1' * Rot(p1_f)' * F];

assert(norm(Q - Q_expected) < 3e-14)



