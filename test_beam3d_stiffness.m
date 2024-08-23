% test for beam3d stiffness matrix

% check is matrix is symmetric and positive semi-definite

E = 2e11;
nu = 0.3;
G = E / 2 / (1 + nu);
l = 0.5;
a = 0.05; % square cross section
A = a * a;
Iy = a ^ 4 / 12;
Iz = Iy;
% Jx = 2 / 9 * a ^ 3; % polar moment of inertia
Jx = a ^ 4 / 6; % polar moment of inertia

K = beam3d_stiffness(E, G, A, l, Jx, Iy, Iz);
tol = eps(K(1, 1)) * 10;

assert(issymmetric(K))

% Check eigen-frequencies
d = eig(K);
% first 6 should be almost zero 
assert(norm(d(1:6)) < tol)
% next 6 should be positive and nonzero
assert( all(d(7:end) > tol) )

%% Force at the end in clambed beam

Fy = -1000;
% Nonzero solution from ANSYS Beam4 model
uy = -0.4e-3; 
rotz = -0.12e-2; 

K_stat = K(7:end, 7:end); % apply boundary conds - clamp at the first node
F = zeros(6, 1);
F(2) = Fy;

x = K_stat\F;

x_expected = [0; uy; 0; 0; 0; rotz];
assert( norm(x - x_expected) < 1e-15 )