% Script to verify how mass matrix looks like

%% Test with one body
bodies.mass = 12.1;
bodies.Ic = rand(3, 3);

sys = sys_generate(bodies, []);

expected = blkdiag(eye(3) * bodies.mass, bodies.Ic);

result = mbd_mass(sys, []);

assert(norm(result - expected) < 1e-15)

%% Three bodies

bodies.mass = 1;
bodies.Ic = eye(3) .* 1.2;
bodies(2).mass = 2;
bodies(2).Ic = eye(3) .* 3.7;
bodies(3).mass = 3;
bodies(3).Ic = eye(3) .* 7.9;

sys = sys_generate(bodies, []);

expected = diag([1, 1, 1, 1.2, 1.2, 1.2, 2, 2, 2, 3.7, 3.7, 3.7, ...
    3, 3, 3, 7.9, 7.9, 7.9]);

result = mbd_mass(sys, []);

assert(norm(result - expected) < 1e-15)
