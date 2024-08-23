%% Test for ffr_grav_force - verify if integral provide the force.

f = h_fbody_integrate(1, 4);

g = [0; 0; -9.81];

p = rand(4, 1);
p = p ./ norm(p);

qf = rand(size(f.Kff, 1), 1) * 0.173174;

Qg = ffr_grav_force(f, g, p, qf);

qn = f.q0 + f.B2 * qf;
Qg_expected = zeros(size(Qg));
% Loop over elements!
for ii = 1 : f.n_elem
    mi = @(xi, eta, zeta) integral_Qg(xi, eta, zeta, f.B2, qn, p, f.elem_len, f.elem_idx(:, ii), g);
    Qg_expected = Qg_expected + volumetric_integral(mi, f.mat, f.sec, f.elem_len, 2, 1, 1);
end

diffQg = Qg_expected - Qg;
assert(norm(diffQg) < 5e-14)

function Qg_int = integral_Qg(xi, eta, zeta, B2, qn, p, elem_len, elem_idx, grav)
    S = zeros(3, size(B2, 1));
    S_elem = beam3d_shape_fun(xi, eta, zeta, elem_len);
    S(:, elem_idx) = S_elem;
    u = S * qn;
    A = Rot(p);
    
    L = [eye(3), -A * sksym(u), A * S * B2];
    Qg_int = L' * grav;
end