% test if beam3d grav procedure is proper

rho = 7801;
l = 0.5;
a = 0.05; % square cross section
A = a * a;
Iy = a ^ 4 / 12;
Iz = Iy;
% Jx = 2 / 9 * a ^ 3; % polar moment of inertia
Jx = a ^ 4 / 6; % polar moment of inertia

grav = [1.2, -0.7, 0.132]';

%% Compare results with triple integral using Gaussian quadrature

g_int = @(xi, eta, zeta) grav_integral(xi, eta, zeta, l, grav);

half_eta = a / 2 / l;
% Tested that ranks 4 and 2 are required for good results
Fg_int = tripleGaussQuad(g_int, 0, 1, -half_eta, half_eta, -half_eta, half_eta, ...
    4, 2, 2);
Fg_expected = Fg_int * rho * l ^ 3; % l^3 due to coordinate transformation

Fg = beam3d_grav(rho, A, l, grav);

diffFg = Fg - Fg_expected;
error = norm(diffFg);
assert( error < 1e-14 )

function grav_int = grav_integral(xi, eta, zeta, l, grav) 
    S = beam3d_shape_fun(xi, eta, zeta, l);
    grav_int = S'*grav;
end