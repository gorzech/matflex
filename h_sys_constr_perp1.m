function sys = h_sys_constr_perp1(sys, q0, body_i_id, body_j_id, v_i, v_j, location)
%H_SYS_CONSTR_PERP1 Helper to generate data for perpendicular constraint of
%type 1. 
v_i = normunit(v_i);
v_j = normunit(v_j);
assert(dot(v_i, v_j) < 1e-5, 'Vectors seems not to be perpendicular!');

p.nc = 1;

p.body_i_id = body_i_id;
p.body_j_id = body_j_id;

[b, isrigid] = h_get_body(sys, body_i_id);
p_b = q0(b.q_idx(4 : 7));
p.v_i = Rot(p_b)' * v_i;
p.psi_i = [];
if ~isrigid
    [~, psi] = u_get_node_id(b, q0, location);
    assert(~isempty(psi), 'Error when looking for node on body %d', body_i_id)
    p.psi_i = psi;
end

[b, isrigid] = h_get_body(sys, body_j_id);
p_b = q0(b.q_idx(4 : 7));
p.v_j = Rot(p_b)' * v_j;
p.psi_j = [];
if ~isrigid
    [~, psi] = u_get_node_id(b, q0, location);
    assert(~isempty(psi), 'Error when looking for node on body %d', body_j_id)
    p.psi_j = psi;
end

sys.nconstr = sys.nconstr + p.nc;
sys.constr.perp1 = [sys.constr.perp1, p];

end

