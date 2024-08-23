% Test if basic call to mbd_quadratic_inertia works as expected

%% Use two bodies

Ic1 = rand(3, 3);
Ic2 = diag([1.2, 0.13, 78e-3]);
bodies.mass = 3.2;
bodies.Ic = Ic1;
bodies(2).mass = 2;
bodies(2).Ic = Ic2;

om1 = rand(3, 1);
om2 = [0; 1e-2; 0.1234112];
h = [0; 0; 0; om1; 1; 0; 0; om2];

sys = sys_generate(bodies, []);

expected_b = [zeros(3,1)
    sksym(om1) * Ic1 * om1
    zeros(3, 1)
    sksym(om2) * Ic2 * om2];

result_b = mbd_quadratic_inertia(sys, 0, h);

assert(norm(expected_b - result_b) < 1e-15)