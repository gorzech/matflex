% Perform modal analysis for system of few beams. To verify if system mass
% and stiffness matrices looks appropriate

% Define material mat and section properties sec
mat.rho = 7801;
l = 0.5;
a = 0.05; % square cross section
sec.A = a * a;
sec.Iy = a ^ 4 / 12;
sec.Iz = sec.Iy;
% Jx = 2 / 9 * a ^ 3; % polar moment of inertia
sec.Jx = a ^ 4 / 6; % polar moment of inertia

mat.E = 2e11;
nu = 0.3;
mat.G = mat.E / 2 / (1 + nu);

n_elem = 16;
elem_len = l / n_elem;
[~, elem_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);
K = fem_stiffness(elem_idx, mat, sec, elem_len);

%% Free-free modes

% Solution from Ansys for 16 beam4 elements (frequencies 7 to 16)
reference_freq = [1020.1
    1020.1
    2747.3
    2747.3
    3145.2
    5071.5
    5232.6
    5232.6
    6320.8
    8362.0];

d = eig(full(K), full(M));

tol = eps(K(1, 1) / M(1, 1)) * 10;
% first 6 should be almost zero 
assert(norm(d(1:6)) < tol)

freq = sqrt(d(7:16)) / 2 / pi;
assert(all(abs(reference_freq - freq) < 0.1))

%% Modal analysis without rigid body modes

% boundary conditions
bc_nodes = [1, 5, 17];
bc_dofs = {[1, 2, 4], 3, [5, 6]};

% Ansys solution for similar model with beam4 elements (first ten freq)
reference_freq = [114.68
    197.46
    1023.8
    1259.7
    1570.7
    2477.3
    2532.7
    2799.1
    4727.3
    5362.4];

Kb = fem_bc_split(bc_nodes, bc_dofs, K);
Mb = fem_bc_split(bc_nodes, bc_dofs, M);
d = eig(full(Kb), full(Mb));

freq = sqrt(d(1:10)) / 2 / pi;
assert(all(abs(reference_freq - freq) < 0.1))