% test sys_ode with rigid body free fall

%% Just single body

bodies = struct('mass', 1.17893);
bodies.Ic = eye(3);

g = 9.81;
grav = [0; 0; -g];

sys = sys_generate(bodies, [], 'grav', grav);

y0 = [1; 0; 2
    1; 0; 0; 0
    0; 0; 0
    0; 0; 0];

of = @(t, y) sys_ode(sys, t, y);

[T, Y] = ode45(of, 0:1e-2:1, y0);

tend = T(end);
exp_y_end = y0;
exp_y_end(3) = exp_y_end(3) - 0.5 * g * tend ^ 2;
exp_y_end(10) = exp_y_end(10) - g * tend;

assert(norm(Y(end, :)' - exp_y_end) < 1e-14)

%% Rotation with constant angular velocity

bodies = struct('mass', 1.17893);
bodies.Ic = eye(3);

grav = [0; 0; 0];

sys = sys_generate(bodies, [], 'grav', grav);

of = @(t, y) sys_ode(sys, t, y);

y0 = [1; 0; 2
    1; 0; 0; 0
    0; 0; 0
    0; 0; 0];

for ii = 11 : 13
    y0(ii) = 2 * pi;
    [~, Y] = ode45(of, 0:1e-2:1, y0);
    
    % tend = T(end);
    exp_y_end = y0;
    cmp_idx = [1:3 8:length(y0)];
    
    t_idx = 51; % half of the motion - should rotate by 180 deg
    assert(norm(Y(t_idx, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)

    A_expected = RotAxis(pi, ii - 10); % Rotation by 180 deg
    A_idx = Rot(Y(t_idx, 4:7)'./norm(Y(t_idx, 4:7)'));
    assert(norm(A_idx - A_expected) < 1e-7)

    t_idx = 101; % end
    assert(norm(Y(t_idx, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)
    % As EP are not unique - check if their rotational matrices are the same
    A_expected  = eye(3); % body should perform full rotation
    A_idx = Rot(Y(t_idx, 4:7)'./norm(Y(t_idx, 4:7)'));
    assert(norm(A_idx - A_expected) < 1e-6)
    y0(ii) = 0;
end
%% Add second body, translational and rotational velocities

bodies = struct('mass', 1.17893);
bodies.Ic = eye(3);
bodies(2).mass = 289.9872;
bodies(2).Ic = diag([7.12, 0.23, 0.87]);

g = 9.81;
grav = [0; 0; -g];

sys = sys_generate(bodies, [], 'grav', grav);

p2 = rand(4, 1);
p2 = p2 / norm(p2);

y0 = [1; 0; 2 % body 1 r
    1; 0; 0; 0 % body 1 p
    0; 0; 0 % body 2 r
    p2 % body 2 p
    0; 0; 0 % body 1 v
    0; 0; 0 % body 1 om
    0; 10; 0 % body 2 v
    0.1; 0; 0.3 % body 2 om
    ];

of = @(t, y) sys_ode(sys, t, y);

[T, Y] = ode45(of, 0:1e-2:1, y0);

tend = T(end);
exp_y_end = y0;
cmp_idx = [1:10 15:23]; % what we can compare directly against each other
exp_y_end(3) = exp_y_end(3) - 0.5 * g * tend ^ 2;
exp_y_end(9) = exp_y_end(9) + 10 * tend; % due to initial velocity
exp_y_end(10) = exp_y_end(10) - 0.5 * g * tend ^ 2;
exp_y_end(17) = exp_y_end(17) - g * tend;
exp_y_end(23) = exp_y_end(23) - g * tend;

assert(norm(Y(end, cmp_idx)' - exp_y_end(cmp_idx)) < 1e-14)

% Check if EP are normalized and angular kinetic energy preserved for body
% 2
p2 = Y(end, 11:14)';
assert(abs(p2' * p2 - 1) < 1e-9)

om2 = Y(end, 24:26)';
Ea_end = 0.5 * om2' * bodies(2).Ic * om2;
om2_0 = Y(1, 24:26)';
Ea_0 = 0.5 * om2_0' * bodies(2).Ic * om2_0;

assert(abs(Ea_end - Ea_0) < 2e-9)
