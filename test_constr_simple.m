% Test for cnstr_simple

% Most basic - one rigid body - check fields, etc.
rng(3)
body = h_body(2, [0.1, 1e-3, 1e-3]);

sys = sys_generate(body, []);
q0 = [0; 0; 0; 1; 0; 0; 0];

sys = h_sys_constr_simple(sys, q0, 1);

assert(sys.nconstr == 6, "Number of constr set incorrectly")
assert(length(sys.constr.simple) == 1)
assert(sys.constr.simple.nc == 6)
assert(all(size(sys.constr.simple.q0) == [6, 1]))

[c, c_q, c_p, gam] = constr_simple(sys, sys.constr.simple, q0, zeros(6, 1));

assert(length(c) == 6)
assert(norm(c) < 1e-16);
assert(norm(c_p) < 1e-16);
assert(norm(gam) < 1e-16);
assert(norm(c_q - diag([1, 1, 1, 0.5, 0.5, 0.5])) < 1e-16);

%% Set constraints on 2 bodies

p2 = rand(4, 1);
p2 = p2 ./ norm(p2);

fbody = h_fbody(1, 2);

q0 = [0.1; 0.3; 0.2; p2; zeros(3, 1); -p2; zeros(fbody.nflex, 1)];

sys = sys_generate(body, fbody); % start fresh

sys = h_sys_constr_simple(sys, q0, 1);
sys = h_sys_constr_simple(sys, q0, 2);

assert(sys.constr.simple(1).nc == 6)
assert(sys.constr.simple(2).nc == 6)

[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(norm(c) < 1e-16);

q0(2) = 0;
[c, ~, ~, ~] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(abs(norm(c) - 0.3) < 1e-16);

%% Now verify derivatives of the constraint

p2 = rand(4, 1);
p2 = p2 ./ norm(p2);

fbody = h_fbody(1, 2);

q0 = [0.1; 0.3; 0.2; p2; zeros(3, 1); -p2; zeros(fbody.nflex, 1)];

sys = sys_generate(body, fbody); % start fresh

sys = h_sys_constr_simple(sys, q0, 1);
sys = h_sys_constr_simple(sys, q0, 2);

assert(sys.nconstr == 12)

h = rand(sys.nh, 1) * 0.7623;
q0 = q0 + 1e-6;
% normalize quaternions
q0(4 : 7) = q0(4 : 7) ./ norm(q0(4 : 7));
q0(11 : 14) = q0(11 : 14) ./ norm(q0(11 : 14));

[~, c_q, c_p, g] = mbd_constr(sys, 0, q0, h);
[expected_c_h, expected_c_p, expected_g] = th_approx_constr_derivatives(sys, @mbd_constr, q0, h);

assert(norm(c_q - expected_c_h) < 6e-7)
assert(norm(c_p - expected_c_p) < 1e-14)
assert(norm(g - expected_g) < 1e-14)

%% Verify when not all rigit parts are constraints 

p2 = rand(4, 1);
p2 = p2 ./ norm(p2);

fbody = h_fbody(1, 2);

q0 = [0.1; 0.3; 0.2; p2; zeros(3, 1); -p2; zeros(fbody.nflex, 1)];

sys = sys_generate(body, fbody); % start fresh

sys = h_sys_constr_simple(sys, q0, 1, [1, 2]);
sys = h_sys_constr_simple(sys, q0, 2, [2, 3]);

assert(sys.nconstr == 10)

[c] = mbd_constr(sys, 0, q0, zeros(sys.nh, 1));
assert(length(c) == 10)
assert(norm(c) < 1e-16);

h = rand(sys.nh, 1) * 0.7623;
q0 = q0 + 1e-6;
% normalize quaternions
q0(4 : 7) = q0(4 : 7) ./ norm(q0(4 : 7));
q0(11 : 14) = q0(11 : 14) ./ norm(q0(11 : 14));

[~, c_q, c_p, g] = mbd_constr(sys, 0, q0, h);
[expected_c_h, expected_c_p, expected_g] = th_approx_constr_derivatives(sys, @mbd_constr, q0, h);

assert(norm(c_q - expected_c_h) < 6e-7)
assert(norm(c_p - expected_c_p) < 1e-14)
assert(norm(g - expected_g) < 1e-14)

