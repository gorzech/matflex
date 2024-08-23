function fbody = h_fbody_integrate(body_len, n_elem, varargin)
%H_FBODY Heleper to generate flexible body - suitable for integrate code
%parts

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'body_len',validScalarPosNum);
addOptional(p, 'n_elem', 2, validScalarPosNum);
addParameter(p, 'sec', h_section());
addParameter(p, 'mat', h_material());
addParameter(p, 'rc_nodes', 1);
addParameter(p, 'rc_dofs', { 1 : 6 });

parse(p, body_len, n_elem, varargin{:});

n_elem = p.Results.n_elem;
elem_len = body_len / n_elem;
[q0, elem_idx] = straight_beam3d(elem_len, n_elem);

M = fem_mass(elem_idx, p.Results.mat, p.Results.sec, elem_len);
K = fem_stiffness(elem_idx, p.Results.mat, p.Results.sec, elem_len);

%   fbody contains (at lest)
fbody.q0 = q0;
fbody.elem_idx = elem_idx;
fbody.mat = p.Results.mat;
fbody.sec = p.Results.sec;
fbody.elem_len = elem_len;
fbody.n_elem = n_elem;
fbody.mass = fbody.mat.rho * fbody.sec.A * elem_len * n_elem;

%   B2 - reference conditions
B2 = speye(size(M));
[~, constr_dof] = fem_boundary_conditions(p.Results.rc_nodes, p.Results.rc_dofs, length(q0));
B2(:, constr_dof) = [];
fbody.B2 = B2;

fbody.m_ff = B2' * M * B2;

fbody.nflex = length(q0) - length(constr_dof);
fbody.nh = 6 + fbody.nflex;
fbody.nq = fbody.nh + 1;

fbody.Kff = B2' * K * B2;

end

