function B2 = h_set_B2(nqf, rc_nodes, rc_dofs)
%H_SET_B2 Set B2 reference conditions
%   B2 - reference conditions
B2 = speye(nqf);
[~, constr_dof] = fem_boundary_conditions(rc_nodes, rc_dofs, nqf);
B2(:, constr_dof) = [];
end

