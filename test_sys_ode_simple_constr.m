% System with simple constr should not move

rng(3);
body = h_body(2, [0.1, 1e-3, 1e-3]);
fbody = h_fbody(1, 2);

p2 = rand(4, 1);
p2 = p2 ./ norm(p2);

q0 = [0.1; 0.3; 0.2; p2; zeros(3, 1); [0.99; 0; 0.14107; 0]; zeros(fbody.nflex, 1)];

%% Try to integrate with one constraint

sys = sys_generate(body, fbody, 'grav', [0; 10; 0], 'usequadratic', false); 
sys = h_sys_constr_simple(sys, q0, 1);

y0 = [q0; zeros(sys.nh, 1)];

of = @(t, y) sys_ode(sys, t, y);
opts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9);
[T, Y] = ode15s(of, 0:1e-2:1, y0, opts);

% Change only for flexible body
tend = T(end);
y_exp = y0;
y_exp(9) = y_exp(9) + 0.5 * 10 * tend ^ 2;
y_exp(sys.nq + 8) = y_exp(sys.nq + 8) + 10 * tend;

assert(norm(Y(end, :)' - y_exp) < 5e-9)

%% Now for flexible
sys = sys_generate(body, fbody, 'grav', [10; 0; 0], 'usequadratic', false); 
sys = h_sys_constr_simple(sys, q0, 1);
sys = h_sys_constr_simple(sys, q0, 2);

y0 = [q0; zeros(sys.nh, 1)];

of = @(t, y) sys_ode(sys, t, y);
opts = odeset('AbsTol', 1e-3, 'RelTol', 1e-4);
% tic
[T, Y] = ode15s(of, 0:1e-2:1, y0, opts);
% toc
% Change only for flexible body
tend = T(end);
% Position and velocity should be constant for the frame
test_idx = [1:14, sys.nq + (1:12)];
assert(norm(Y(end, test_idx)' - y0(test_idx)) < 1e-7)
