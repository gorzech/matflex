% Perform modal analysis for system of few beams. To verify if system mass
% and stiffness matrices looks appropriate

% Define material mat and section properties sec
mat.rho = 7801;
l = 2;
% Thin-walled cross section
sec.A = 0.006144;
sec.Iy = 3.7814e-05;
sec.Iz = 3.7814e-05;
% This test verify results when Jx != Iy = Iz
sec.Jx = 4.9836032e-05; % sec.Iy + sec.Iz;

mat.E = 2e11;
nu = 0.3;
mat.G = mat.E / 2 / (1 + nu);

n_elem = 40;
elem_len = l / n_elem;
[~, elem_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);
K = fem_stiffness(elem_idx, mat, sec, elem_len);

% %% Free-free modes
% 
% % Solution from Ansys for 16 beam4 elements (frequencies 7 to 16)
% reference_freq = [1020.1
%     1020.1
%     2747.3
%     2747.3
%     3145.2
%     5071.5
%     5232.6
%     5232.6
%     6320.8
%     8362.0];
% 
% d = eig(full(K), full(M));
% 
% tol = eps(K(1, 1) / M(1, 1)) * 10;
% % first 6 should be almost zero 
% assert(norm(d(1:6)) < tol)
% 
% freq = sqrt(d(7:16)) / 2 / pi;
% assert(all(abs(reference_freq - freq) < 0.1))

%% Modal analysis without rigid body modes

% boundary conditions
bc_nodes = [1, 41];
bc_dofs = {1:6, 1:6};

% Ansys solution for similar model with beam4 elements (first ten freq)
reference_freq = [ 350.31
     350.31
     637.43
     941.92
     941.92
     1266.2
     1275.9
     1780.2
     1780.2
     1916.2
     2534.3
     2559.6
     2809.5
     2809.5
     3206.9
     3806.3
     3859.1
     3979.2];

Kb = fem_bc_split(bc_nodes, bc_dofs, K);
Mb = fem_bc_split(bc_nodes, bc_dofs, M);
d = eig(full(Kb), full(Mb));

freq = sqrt(d(1:18)) / 2 / pi;
assert(all(abs(reference_freq - freq) < 0.1))