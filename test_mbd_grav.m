% Short test for gravity vector computations

%% 2 bodies

g = [0.1; -9.81; 6];

m1 = 3.2;
m2 = 2;
Ic1 = rand(3, 3);
Ic2 = diag([1.2, 0.13, 78e-3]);
bodies.mass = m1;
bodies.Ic = Ic1;
bodies(2).mass = m2;
bodies(2).Ic = Ic2;

sys = sys_generate(bodies, [], 'grav', g);

expected_g = [m1 .* g
    zeros(3, 1)
    m2 .* g
    zeros(3, 1)];

result_g = mbd_grav(sys, 0);

assert(norm(expected_g - result_g) < 1e-15)