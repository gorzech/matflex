% Simple test for boundary conditions

%% Really simple test
M = reshape(1:36, [6, 6]);
node_id = 1;
dofs = {1};

[Mi, Mb, Mib, ~, P]  = fem_bc_split(node_id, dofs, M);
assert(all(size(Mi) == [5, 5]))
assert(all(Mi(:, 1) == (8:12)'))
assert(Mb == 1)
assert(all(Mib == (2:6)'))
assert(all(P == (1:6)'))

%% Test with dofs removed from one node
K = rand(18);

node_id = 2; 
dofs = {[1 2 3 5]}; % global dofs [7, 8, 9, 11]

[Ki, Kb, Kib, Kbi, P]  = fem_bc_split(node_id, dofs, K);
assert(all(size(Ki) == [14, 14]))
assert(all(size(Kb) == [4, 4]))
assert(all(size(Kib) == [14, 4]))
assert(all(size(Kbi) == [4, 14]))

b_dofs = [7, 8, 9, 11];
i_dofs = [1:6, 10, 12:18];
assert( all(all(Ki == K(i_dofs, i_dofs))) )
assert( all(all(Kb == K(b_dofs, b_dofs))) )
assert( all(all(Kbi == K(b_dofs, i_dofs))) )
assert( all(all(Kib == K(i_dofs, b_dofs))) )
assert( all(all(K(P, P) == [Kb, Kbi; Kib, Ki])) )

%% Test with multiple nodes
K = rand(18);

node_id = 1:3;
dofs = {[1, 3], 1, 3}; % globals [1, 3, 7, 15]

Ki = fem_bc_split(node_id, dofs, K);
assert(all(size(Ki) == [14, 14]))

K([1, 3, 7, 15], :) = [];
K(:, [1, 3, 7, 15]) = [];
assert(all(Ki(:) == K(:)))