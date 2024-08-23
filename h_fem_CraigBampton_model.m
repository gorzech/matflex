function [V, freq] = h_fem_CraigBampton_model(K, M, varargin)
%H_FEM_CRAIGBAMPTON_MODEL Return transformation matrix for the CB method

p = inputParser;
validScalarNonNegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addParameter(p, 'n_modal', 0, validScalarNonNegNum);
addParameter(p, 'cb_nodes', 1);
addParameter(p, 'cb_dofs', { });
addParameter(p, 'two_step', true);
addParameter(p, 'rc_mean_axis', true);
addParameter(p, 'rc_nodes', 1);
addParameter(p, 'rc_dofs', { });
addParameter(p, 'stats', false);

parse(p, varargin{:});

two_step = p.Results.two_step;
mean_axis = p.Results.rc_mean_axis;
% if mean_axis is set to default
if any(strcmp(p.UsingDefaults, 'rc_mean_axis'))
    mean_axis = ~~two_step;
end
if mean_axis
    assert(two_step, "Mean axis reference condition is valid only for two step CB");
end

cb_nodes = p.Results.cb_nodes;
cb_dofs = fill_node_dofs(p.Results.cb_dofs, length(cb_nodes));
    
if two_step
    [V, ~, freq] = fem_CraigBampton_orth(p.Results.n_modal, K, M, ...
        cb_nodes, cb_dofs);
else
    V = fem_CraigBampton(p.Results.n_modal, K, M, cb_nodes, ...
        cb_dofs);
    freq = [];
end

if mean_axis
    tol = sqrt(eps(K(1, 1) / M(1, 1))) * 10;
    rb_modes = sum(freq < tol);
    if p.Results.stats
        fprintf('\n\tYour system has %d RB modes!\n\n', rb_modes);
    end
    if rb_modes ~= 6
        warning('Model might be incorrect as the number of rigid modes is not equal to 6!');
    end
    V = V(:, rb_modes + 1 : end);
else
    rc_nodes = p.Results.rc_nodes;
    rc_dofs = fill_node_dofs(p.Results.rc_dofs, length(rc_nodes));
    [~, constr_dof] = fem_boundary_conditions(rc_nodes, rc_dofs, size(M, 1));
    V(constr_dof, :) = 0; % exclude all those boundary stuff. Should be also interface nodes
end

end

function cb_dofs = fill_node_dofs(curr_dofs, n)
cb_dofs = cell(1, n);
for ii = 1 : n
    if ii > length(curr_dofs)
        cb_dofs{ii} = 1:6;
    elseif isempty(curr_dofs{ii})
        cb_dofs{ii} = 1:6;
    else
        cb_dofs{ii} = curr_dofs{ii};
    end
end
end