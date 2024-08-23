function fbody = h_fbody(body_len, n_elem, varargin)
%H_FBODY Heleper to generate flexible body - current version that does work
%with preintegrated scripts

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'body_len',validScalarPosNum);
addOptional(p, 'n_elem', 2, validScalarPosNum);
addParameter(p, 'sec', h_section());
addParameter(p, 'mat', h_material());
addParameter(p, 'rc_nodes', 1);
addParameter(p, 'rc_dofs', { 1 : 6 });
addParameter(p, 'cb_opts', []);
addParameter(p, 'move_to_cm', false);
addParameter(p, 'lumped_mass', false);
addParameter(p, 'r_alpha', 0);
addParameter(p, 'r_beta', 0);

parse(p, body_len, n_elem, varargin{:});

[qf0, K, M, node_idx] = h_fem_model(body_len, p.Results.n_elem, ...
    p.Results.sec, p.Results.mat, p.Results.lumped_mass);

cb_opts = p.Results.cb_opts;
if isempty(cb_opts)
    B2 = h_set_B2(size(M, 1), p.Results.rc_nodes, ...
        p.Results.rc_dofs);
else
    cb_opts.rc_nodes = p.Results.rc_nodes;
    cb_opts.rc_dofs = p.Results.rc_dofs;
    B2 = h_fem_CraigBampton_model(K, M, cb_opts);
end

fbody = h_fbody2(qf0, M, K, B2, node_idx, p.Results.move_to_cm);

if p.Results.r_alpha > 0 || p.Results.r_beta > 0
    fbody.De = p.Results.r_alpha .* fbody.m_ff + ...
        p.Results.r_beta .* fbody.Kff;
end

end

