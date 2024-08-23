function sys = h_sys_constr_simple(sys, q0, body_id, r_idx_fix)
%H_SYS_CNSTR_SIMPLE Function to generate base for simple constraint that
%fixes coordinate system of the body

if nargin < 4
    r_idx_fix = 1:3;
end
r_idx_fix = r_idx_fix(:)'; % make row vector

% Need to set nc, Bq0, B, h_idx (global)
body = h_get_body(sys, body_id);
s.h_idx = body.h_idx([r_idx_fix, 4 : 6]);
s.nc = length(s.h_idx);

s.q_idx = body.q_idx([r_idx_fix, 5 : 7]);
s.q0 = q0(s.q_idx);
s.G0_2 = 0.5 .* Gep(q0(body.q_idx(4 : 7)));
s.G0_2 = s.G0_2(:, 2 : end)';

sys.nconstr = sys.nconstr + s.nc;
sys.constr.simple = [sys.constr.simple, s];

end

