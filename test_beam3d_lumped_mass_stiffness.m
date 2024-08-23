% test if beam3d mass & stiffness provide proper solution for some basic
% analysis

rho = 7801;
l = 0.5;
a = 0.05; % square cross section
A = a * a;
Iy = a ^ 4 / 12;
Iz = Iy;
% Jx = 2 / 9 * a ^ 3; % polar moment of inertia
Jx = a ^ 4 / 6; % polar moment of inertia

E = 2e11;
nu = 0.3;
G = E / 2 / (1 + nu);

lumped = true;
M = beam3d_mass(rho, A, l, Jx, Iy, Iz, lumped);
K = beam3d_stiffness(E, G, A, l, Jx, Iy, Iz);

%% Modal analysis for free-free beam model

% modal frequencies (6) from ANSYS Beam4 model with one element
% First 5 are zero
% Rest are undetermined?
reference_freq = [3223.4];

d = eig(K, M);

dr = sort(d(isfinite(d)));

tol = eps(K(1, 1) / M(1, 1)) * 10;
% first few should be almost zero 
assert(norm(dr(1:end - 1)) < tol)

freq = sqrt(dr(end)) / 2 / pi;
assert(all(abs(reference_freq - freq) < 0.1))