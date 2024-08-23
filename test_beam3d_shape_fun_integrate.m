%% Test for beam3d_shape_fun_integrate - verify if integral is correct

mat = h_material();
sec = h_section();

L = 0.63263;

S_int = beam3d_shape_fun_integrate(mat.rho, sec.A, L);

mi = @(xi, eta, zeta) integral_S(xi, eta, zeta, L);
S_expected = volumetric_integral(mi, mat, sec, L, 4, 2, 2);

diffS = S_expected - S_int;
assert(norm(diffS) < 1e-14)

function S = integral_S(xi, eta, zeta, elem_len)
    S = beam3d_shape_fun(xi, eta, zeta, elem_len);
end