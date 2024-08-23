% Try to make basic test for normalized_eig. An easy approach - based on
% beam3d matrices

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

K = beam3d_stiffness(E, G, A, l, Jx, Iy, Iz);

%% Test with full matrices

M = beam3d_mass(rho, A, l, Jx, Iy, Iz);

[Ve, De] = eig(K, M, 'vector');

[V, F] = normalized_eig(K, M, 11);

assert(norm(Ve(:, 1:11) - V) < 1e-13)
assert(norm( F - sqrt(abs(De(1 : 11))) / 2 / pi ) < 1e-16)

%% Second test with lumped mass matrix

M_lumped = beam3d_mass(rho, A, l, Jx, Iy, Iz, true);

reference_freq = [3223.4];

% From other test I assume it returns at most 5 values
% in any case all except the last one are close to zero
[V, F] = normalized_eig(K, M_lumped, 7);

tol = sqrt(eps(K(1, 1) / M_lumped(1, 1))) / 10;
% first few should be almost zero 
assert(norm(F(1 : end - 1)) < tol)

assert(all(abs(reference_freq - F(end)) < 0.1))

% Are they M-normalized?
assert(norm(V' * M_lumped * V - eye(size(V, 2))) < 0.002)

[Ve, De] = eig(K, M_lumped, 'vector');
idx_freq = find(De > 1000 & De < Inf); % there is one finite positive value

% Normalize corresponding vectors to one
Ve_norm = Ve(:, idx_freq) ./ norm(Ve(:, idx_freq));
V_norm = V(:, end) ./ norm(V(:, end));

% Now they should be the same vector
assert(norm(Ve_norm - V_norm) < 1e-15 || norm(Ve_norm > V_norm) < 1e-15)