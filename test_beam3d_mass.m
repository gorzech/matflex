% test if beam3d mass procedure is proper

% Test is mass matrix is symmetric and seems to be positive definite
% If this fails does not proceed further
rho = 7801;
l = 0.5;
a = 0.05; % square cross section
A = a * a;
Iy = a ^ 4 / 12;
Iz = Iy;
% Jx = 2 / 9 * a ^ 3; % polar moment of inertia
Jx = a ^ 4 / 6; % polar moment of inertia

M = beam3d_mass(rho, A, l, Jx, Iy, Iz);

[R, flag] = chol(M);
assert(flag == 0) % if flag is 0 matrix is symmetric positive definite

%% Now compare results with triple integral using Gaussian quadrature

m_int = @(xi, eta, zeta) mass_integral(xi, eta, zeta, l);

half_eta = a / 2 / l;
% Tested that ranks 4 and 2 are required for good results
M_int = tripleGaussQuad(m_int, 0, 1, -half_eta, half_eta, -half_eta, half_eta, ...
    4, 2, 2);
M_result = M_int * rho * l ^ 3; % l^3 due to coordinate transformation

diffM = M - M_result;
error = norm(diffM);
assert( error < 1e-14 )

function m_int = mass_integral(xi, eta, zeta, l) 
    S = beam3d_shape_fun(xi, eta, zeta, l);
    m_int = S'*S;
end