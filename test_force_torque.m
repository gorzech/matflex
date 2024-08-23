%% Rigid body

rng(98);
r1 = [0; 1; -1];
p1 = normunit(rand(4, 1));

sys = sys_generate(h_body(1, eye(3)), [], 'usequadratic', false);

q0 = [r1; p1];

torque_point = [0; 1; -1];
s0 = r1 + [0; 1; -1];

torque = @(t) [3; 2 + 2*t; t^2 - 1];

sys = h_sys_torque(sys, torque, s0, q0, 1);

assert(length(sys.force.torque) == 1);

Q = mbd_force(sys, 0.1, q0, zeros(sys.nh, 1));
assert(length(Q) == 6);

L = torque(0.1);
Q_expected = [zeros(3, 1);
    Rot(p1)' * L];

assert(norm(Q - Q_expected) < 1e-14)

%% Check if torque can support body

m = 1.25;
l = 0.7;
g = 9.81;

% horizontal beam under gravity in y 
% equivalent torque
torque = @(t) [0; 0; m * g * l/2];

ground = h_body(1, eye(3));
horizontal_beam = h_body(m, (m * l ^ 2) * [1 / 1200; 1 / 12; 1 / 12]);

sys = sys_generate([ground, horizontal_beam], [], 'grav', [0; -g; 0]);

p0 = [1; 2; 3; 0];
p0 = p0 ./ norm(p0);
q0 = [zeros(3, 1); p0; l / 2; 0; 0; p0];

sys = h_sys_constr_simple(sys, q0, 1); % Fix the ground
s_rot = [0; 0; 0];
rot_axis = [0; 0; 1];
sys = h_sys_constr_point(sys, s_rot, q0, 1, 2);
sys = h_sys_constr_perp1(sys, q0, 1, 2, rot_axis, [1; 0; 0], s_rot);
sys = h_sys_constr_perp1(sys, q0, 1, 2, rot_axis, [0; 1; 0], s_rot);
% Add equivalent torque
sys = h_sys_torque(sys, torque, s_rot, q0, 2);

of = @(t, y) sys_ode(sys, t, y);

y0 = [q0; zeros(sys.nh, 1)];
[T, Y] = ode15s(of, [0, 1], y0);

y_end_err = norm(Y(end, :)' - y0);
assert(y_end_err < 1e-14);

%% Check if torque can support flexible body? When it is reasonably stiff this should work

m = 1.25;
l = 0.7;
g = 9.81;

% horizontal beam under gravity in y 

ground = h_body(1, eye(3));
n_elem = 4;
cb_opts.cb_nodes = [1, n_elem + 1];
% cb_opts.cb_dofs = {[], 1};
mat = h_material(1000, 2e11);
a = sqrt(m / mat.rho / l);
sec = h_section(a);
horizontal_beam = h_fbody(l, n_elem, 'mat', mat, 'cb_opts', cb_opts, ...
    'sec', sec, 'move_to_cm', true);

horizontal_beam.De = 3 * horizontal_beam.m_ff + 1e-4 * horizontal_beam.Kff;

sys = sys_generate(ground, horizontal_beam, 'grav', [0; -g; 0]);

% p0 = [1; 2; 3; 0];
% p0 = p0 ./ norm(p0);
% p0 = q_axis(pi, 3);
p0 = [1; 0; 0; 0];

s_rot = [0; 0; 0];

r_flex_0 = s_rot - Rot(p0) * [-l / 2; 0; 0];
q0 = [zeros(3, 1); p0; r_flex_0; p0; zeros(horizontal_beam.nflex, 1)];

% equivalent torque in oposite direction
torque = @(t) sksym(r_flex_0) * [0; m * g; 0];

sys = h_sys_constr_simple(sys, q0, 1); % Fix the ground

rot_axis = [0; 0; 1];
sys = h_sys_constr_point(sys, s_rot, q0, 1, 2);
sys = h_sys_constr_perp1(sys, q0, 1, 2, rot_axis, [1; 0; 0], s_rot);
sys = h_sys_constr_perp1(sys, q0, 1, 2, rot_axis, [0; 1; 0], s_rot);
% Add equivalent torque
sys = h_sys_torque(sys, torque, s_rot, q0, 2);

of = @(t, y) sys_ode(sys, t, y);

y0 = [q0; zeros(sys.nh, 1)];
opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-3);
[T, Y] = ode15s(of, [0, 1e-1], y0);

q_end_err = norm(Y(end, 1 : sys.nq)' - q0);
assert(q_end_err < 2e-6);