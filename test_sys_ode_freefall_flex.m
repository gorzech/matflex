%% Freefall of the flexible body
rng(15)
fbody = h_fbody(1, 4);

g = 9.81;
grav = [0; 0; -g];
sys = sys_generate([], fbody, 'grav', grav);

p2 = rand(4, 1);
p2 = p2 / norm(p2);

y0 = [1; 0; 2
    p2
    zeros(fbody.nflex, 1)
    zeros(fbody.nh, 1)
    ];

of = @(t, y) sys_ode(sys, t, y);

opts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9); %, 'OutputFcn', @odeplot, 'OutputSel', 3);
[T, Y] = ode15s(of, [0, 1], y0, opts);

tend = T(end);
exp_y_end = y0;
exp_y_end(3) = exp_y_end(3) - 0.5 * g * tend ^ 2;
exp_y_end(fbody.nq + 3) = exp_y_end(fbody.nq + 3) - g * tend;

% Error on position
assert(norm(Y(end, 1:fbody.nq)' - exp_y_end(1:fbody.nq)) < 1e-8)
% Error on velocity
assert(norm(Y(end, fbody.nq+1:end)' - exp_y_end(fbody.nq+1:end)) < 1e-13)

%% Same as above but for the integrated version of it
rng(30)
fbody = h_fbody_integrate(1, 4);

g = 9.81;
grav = [0; 0; -g];
sys = sys_generate([], fbody, 'grav', grav);

p2 = rand(4, 1);
p2 = p2 / norm(p2);

y0 = [1; 0; 2
    p2
    zeros(fbody.nflex, 1)
    zeros(fbody.nh, 1)
    ];

of = @(t, y) sys_ode(sys, t, y);

opts = odeset('AbsTol', 1e-7, 'RelTol', 1e-7); %, 'OutputFcn', @odeplot, 'OutputSel', 3);
[T, Y] = ode15s(of, [0, 1], y0, opts);

tend = T(end);
exp_y_end = y0;
exp_y_end(3) = exp_y_end(3) - 0.5 * g * tend ^ 2;
exp_y_end(fbody.nq + 3) = exp_y_end(fbody.nq + 3) - g * tend;

% Error on position
assert(norm(Y(end, 1:fbody.nq)' - exp_y_end(1:fbody.nq)) < 5e-7)
% Error on velocity
assert(norm(Y(end, fbody.nq+1:end)' - exp_y_end(fbody.nq+1:end)) < 1e-14)

%% Rotation with constant angular velocity
% Works only when no quadratic velocity vector is included

fbody = h_fbody(1, 4);

sys = sys_generate([], fbody, 'usequadratic', false);

of = @(t, y) sys_ode(sys, t, y);

y0 = [2; 1; 1
    1; 0; 0; 0
    zeros(fbody.nflex, 1)
    zeros(fbody.nh, 1)
    ];

y0(fbody.nq + 5) = 2 * pi; % set om_y to 2 pi

opts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9);
[~, Y] = ode15s(of, [0, 0.5, 1], y0, opts);

exp_y_end = y0;
cmp_idx = [1:3 8:length(y0)];

t_idx = 2; % half of the motion - should rotate by 180 deg
assert(norm(Y(t_idx, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)

A_expected = RotAxis(pi, 2); % Rotation by 180 deg aroung y axis
A_idx = Rot(Y(t_idx, 4:7)'./norm(Y(t_idx, 4:7)'));
assert(norm(A_idx - A_expected) < 1e-7)

t_idx = 3; % end
assert(norm(Y(t_idx, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)
% As EP are not unique - check if their rotational matrices are the same
A_expected  = eye(3); % body should perform full rotation
A_idx = Rot(Y(t_idx, 4:7)'./norm(Y(t_idx, 4:7)'));
assert(norm(A_idx - A_expected) < 1e-7)

%% Rotation with constant angular velocity - integrated version

fbody = h_fbody_integrate(1, 4);

sys = sys_generate([], fbody, 'usequadratic', false);

of = @(t, y) sys_ode(sys, t, y);

y0 = [2; 1; 1
    1; 0; 0; 0
    zeros(fbody.nflex, 1)
    zeros(fbody.nh, 1)
    ];

y0(fbody.nq + 5) = 2 * pi; % set om_y to 2 pi

opts = odeset('AbsTol', 1e-7, 'RelTol', 1e-7);
[~, Y] = ode15s(of, [0, 0.5, 1], y0, opts);

exp_y_end = y0;
cmp_idx = [1:3 8:length(y0)];

t_idx = 2; % half of the motion - should rotate by 180 deg
assert(norm(Y(t_idx, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)

A_expected = RotAxis(pi, 2); % Rotation by 180 deg aroung y axis
A_idx = Rot(Y(t_idx, 4:7)'./norm(Y(t_idx, 4:7)'));
assert(norm(A_idx - A_expected) < 1e-6)

t_idx = 3; % end
assert(norm(Y(t_idx, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)
% As EP are not unique - check if their rotational matrices are the same
A_expected  = eye(3); % body should perform full rotation
A_idx = Rot(Y(t_idx, 4:7)'./norm(Y(t_idx, 4:7)'));
assert(norm(A_idx - A_expected) < 5e-6)

%% Freefall of the CB flexible body
rng(168)
cb_opts.cb_nodes = [1, 5];
fbody = h_fbody(1, 4, 'cb_opts', cb_opts, 'move_to_cm', true);

assert(fbody.nflex == 6)

g = 9.81;
grav = [0; 0; -g];
sys = sys_generate([], fbody, 'grav', grav);

p2 = rand(4, 1);
p2 = p2 / norm(p2);

y0 = [1; 0; 2
    p2
    zeros(fbody.nflex, 1)
    zeros(fbody.nh, 1)
    ];

of = @(t, y) sys_ode(sys, t, y);

opts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9); %, 'OutputFcn', @odeplot, 'OutputSel', 3);
[T, Y] = ode15s(of, [0, 1], y0, opts);

tend = T(end);
exp_y_end = y0;
exp_y_end(3) = exp_y_end(3) - 0.5 * g * tend ^ 2;
exp_y_end(fbody.nq + 3) = exp_y_end(fbody.nq + 3) - g * tend;

% Error on position
assert(norm(Y(end, 1:fbody.nq)' - exp_y_end(1:fbody.nq)) < 1e-8)
% Error on velocity
assert(norm(Y(end, fbody.nq+1:end)' - exp_y_end(fbody.nq+1:end)) < 1e-13)