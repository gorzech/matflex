% Some basic tests for Craig-Bampton system reduction method

% copy paste from system mass stiffness
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
[qf, elem_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);
K = fem_stiffness(elem_idx, mat, sec, elem_len);

% Solution from Ansys for 16 beam4 elements (frequencies 7 to 16)
reference_freq_free = [1020.1
    1020.1
    2747.3
    2747.3
    3145.2
    5071.5
    5232.6
    5232.6
    6320.8
    8362.0];

%% Free-free modes of transformed system - no actual reduction

% boundary conditions for Craig-Bampton
bc_nodes = [1, 17];
bc_dofs = {1:6, 1:6};

[V, Kv] = fem_CraigBampton_orth(15 * 6, K, M, bc_nodes, bc_dofs);

assert(norm(V' * M * V - eye(size(M))) < 1e-12) % Check if mass matrix is almost unity
% coordinates are transformed but the system should remain almost the same

d = sort(eig(Kv));

tol = eps(K(1, 1) / M(1, 1)) * 100; 
% first 6 should be almost zero 
assert(norm(d(1:6)) < tol)

% But highest system frequencies holds
freq = sqrt(d(7:16)) / 2 / pi;
assert(all(abs(reference_freq_free - freq) < 0.1))

%% Change interface dofs and recheck

% boundary conditions for Craig-Bampton
bc_nodes = [4, 8];
bc_dofs = {1:3, 4:6};

[V, Kv] = fem_CraigBampton_orth(16 * 6, K, M, bc_nodes, bc_dofs);

% coordinates are transformed but the system should remain almost the same
d = sort(eig(Kv));

tol = eps(K(1, 1) / M(1, 1)) * 100; % quite cruel approx is required here
% first 6 should be almost zero 
assert(norm(d(1:6)) < tol)

% But highest system frequencies holds
freq = sqrt(d(7:16)) / 2 / pi;
assert(all(abs(reference_freq_free - freq) < 0.1))

%% Some actual reduction with clamped-clamped interface

% boundary conditions for Craig-Bampton
bc_nodes = [1, 17];
bc_dofs = {1:6, 1:6};

n_modes = 12;
[V, Kv] = fem_CraigBampton_orth(n_modes, K, M, bc_nodes, bc_dofs);
assert(all(size(V) == [17 * 6, 12 + n_modes]))
assert(all(size(Kv) == [12 + n_modes, 12 + n_modes]))

d = sort(eig(Kv));

tol = eps(K(1, 1) / M(1, 1)) * 100; % quite cruel approx is required here

% first 6 should be almost zero 
assert(norm(d(1:6)) < tol)

% But highest system frequencies holds
freq = sqrt(d(7:16)) / 2 / pi;
% Less accuracy - require 10 % accuracy at all modes under consideration
assert(all(abs(reference_freq_free - freq) ./ reference_freq_free < 0.1))