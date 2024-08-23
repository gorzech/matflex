function sys = h_sys_constr_point(sys, location, q0, body_i_id, body_j_id)
%H_SYS_CONSTR_POINT Helper function for adding point (spherical) constraint
%parameters

p.nc = 3;
p.body_i_id = body_i_id;
p.body_j_id = body_j_id;

[b, isrigid] = h_get_body(sys, body_i_id);
if isrigid
    p.v_i = u_global_to_local(b, q0, location);
else
    p.v_i = u_get_node_id(b, q0, location);
    assert(~isempty(p.v_i), 'Error when looking for node on body %d', body_i_id)
end

[b, isrigid] = h_get_body(sys, body_j_id);
if isrigid
    p.v_j = u_global_to_local(b, q0, location);
else
    p.v_j = u_get_node_id(b, q0, location);
    assert(~isempty(p.v_j), 'Error when looking for node on body %d', body_j_id)
end

sys.nconstr = sys.nconstr + p.nc;
sys.constr.point = [sys.constr.point, p];

end