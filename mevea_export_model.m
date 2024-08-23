function mevea_export_model(file_prefix, nodes_loc, V, K_modal, nodal_masses)
% Function to export data in format acceptable by Mevea software
%   Input:
%   file_prefix - file name prefix, 
%   nodes_loc - cartesian position of the nodes (x, y, z) in local
%   coordinates n_nodes x 3 matrix
%   V - transformation matrix
%   K_modal - reduced stiffness matrix
%   nodal_masses - nodal masses (as for lumped mass matrix)

n_nodes = size(nodes_loc, 1);
assert(3 == size(nodes_loc, 2), "node_loc must have 3 columns!")
assert(length(nodal_masses) == n_nodes)
n_modes = size(V, 2);
assert(all(size(K_modal) == [n_modes, n_modes]))

% file initial nodes position
file_N = sprintf('%s_N.dat', file_prefix);
headerC = cell(1, 3);
headerC{1} = n_nodes;

cell_N = [headerC; num2cell(nodes_loc)];
writecell(cell_N, file_N, 'Delimiter', ' ');

% file MO - modes
file_MO = sprintf('%s_MO.dat', file_prefix);

headerC = cell(1, n_modes);
headerC{1} = n_modes;
cell_MO = [headerC; num2cell(V)];

writecell(cell_MO, file_MO, 'Delimiter', ' ');

% file K
file_K = sprintf('%s_K.dat', file_prefix);
K_modal(abs(K_modal) < eps(1e3*K_modal(1,1))) = 0;
writecell(num2cell(K_modal), file_K, 'Delimiter', ' ');

% file NM - nodal mass
file_NM = sprintf('%s_NM.dat', file_prefix);
writecell(num2cell(nodal_masses), file_NM, 'Delimiter', ' ');