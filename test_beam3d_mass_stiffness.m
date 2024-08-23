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

M = beam3d_mass(rho, A, l, Jx, Iy, Iz);
K = beam3d_stiffness(E, G, A, l, Jx, Iy, Iz);

%% Modal analysis for free-free beam model

% modal frequencies (7 to 12) from ANSYS Beam4 model with one element
reference_freq = [1218.3
    1218.3
    3462.5
    3994.2
    3994.2
    5583.2];

d = eig(K, M);

tol = eps(K(1, 1) / M(1, 1)) * 10;
% first 6 should be almost zero 
assert(norm(d(1:6)) < tol)

freq = sqrt(d(7:12)) / 2 / pi;
assert(all(abs(reference_freq - freq) < 0.1))