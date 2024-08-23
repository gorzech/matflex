% Perform modal analysis for system of few beams. To verify if system mass
% and stiffness matrices in ffr looks appropriate

% Define material mat and section properties sec
mat.rho = 7801;
l = 0.5;
a = 0.05; % square cross section
sec.A = a * a;
sec.Iy = a ^ 4 / 12;
sec.Iz = sec.Iy;
sec.Jx = a ^ 4 / 6; % polar moment of inertia
sec.thky = a;
sec.thkz = a;

mat.E = 2e11;
nu = 0.3;
mat.G = mat.E / 2 / (1 + nu);

n_elem = 16;
elem_len = l / n_elem;
[q0, elem_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, mat, sec, elem_len);
K = fem_stiffness(elem_idx, mat, sec, elem_len);

%   fbody contains (at lest)
%   B2 - reference conditions
fbody.q0 = q0;
fbody.elem_idx = elem_idx;
fbody.mat = mat;
fbody.sec = sec;
fbody.elem_len = elem_len;
fbody.n_elem = n_elem;
fbody.mass = mat.rho * sec.A * l;
B2 = speye(size(M));
% Clamped floating frame at first node

% boundary conditions
% Simply supported beam - reasonably similar to free-free
bc_nodes = [1, 17];
bc_dofs = { 1:4, 2:3 };
[~, constr_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, length(q0));

B2(:, constr_dof) = [];
fbody.B2 = B2;

fbody.m_ff = B2' * M * B2;

% test call
p = [1; 0; 0; 0];
qf = zeros(length(q0) - 6, 1);
M_ffr0 = ffr_mass_integrate(fbody, p, qf);
assert(issymmetric(M_ffr0))

K_ff = blkdiag(sparse(6, 6), B2' * K * B2);

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

tol = eps(K(1, 1) / M(1, 1)) * 10;

[~, freq] = normalized_eig(K_ff, M_ffr0, 18);
% first 6 should be almost zero 
assert(norm(freq(1:6)) < tol)

% Skip torsion - 1st bnd, 2nd bnd, axial
ref_modes_id = [1, 2, 3, 4, 6, 7, 8, 10];
ffr_modes_id = 6 + [1, 2, 4, 5, 7, 8, 9, 11];

rel_diff_perc = 100 * abs(reference_freq(ref_modes_id) - freq(ffr_modes_id)) ./ reference_freq(ref_modes_id);
assert(all(rel_diff_perc < 1)) % All bending and axial modes difference smaller than 1%
