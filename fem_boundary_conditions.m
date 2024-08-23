function [free_dofs, constr_dof] = fem_boundary_conditions(bc_nodes, bc_dofs, dofs_count)
%FEM_BOUNDARY_CONDITIONS Generate matrices/vectors to aplly boundary
%conditions on FEM structure
%   Input:
%   bc_nodes - list of nodes where to apply conditions - a vector
%   bc_dofs - cell array of size of node_id with associated dofs to be removed
%   (1 to 6 for our beam)
%   dofs_count - total number of dofs
%   Output:
%   free_dofs - dofs that are here to stay
%   constr_dofs - dofs that are associated with BC - constrained dofs

assert(isvector(bc_nodes), 'node_id must be a vector')
assert(iscell(bc_dofs), 'dofs must be a cell array')
assert(length(bc_dofs) == length(bc_dofs), 'node_id and dofs must have same length')

node_dofs = 6;
constr_dof = zeros(node_dofs * length(bc_nodes), 1); % allocate some space
n_constr_dof = 0;


for ii = 1 : length(bc_nodes)
    assert(bc_nodes(ii) * node_dofs <= dofs_count && bc_nodes(ii) > 0, ...
        'Node %d at %d out of range', bc_nodes(ii), ii)
    assert(all( bc_dofs{ii} > 0 & bc_dofs{ii} <= node_dofs ), ...
        'Some dofs at %d are out of range', ii)
    begin_dof = node_dofs * (bc_nodes(ii) - 1);
    ld = length(bc_dofs{ii});
    constr_dof(n_constr_dof + (1 : ld)) = ...
        begin_dof + bc_dofs{ii};
    n_constr_dof = n_constr_dof + ld;
end
constr_dof = constr_dof(1 : n_constr_dof);

free_dofs = (1:dofs_count)';
free_dofs(constr_dof) = [];

end

