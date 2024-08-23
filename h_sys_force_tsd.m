function sys = h_sys_force_tsd(sys, k, d, l0, point_i, point_j, q0, body_i_id, body_j_id)
%H_SYS_FORCE_TSD Helper to add tsd between two bodies

p.k = k;
p.d = d;
p.l0 = l0;

p.body_i_id = body_i_id;
p.body_j_id = body_j_id;

[b, isrigid] = h_get_body(sys, body_i_id);
if isrigid
    p.v_i = u_global_to_local(b, q0, point_i);
else
    p.v_i = u_get_node_id(b, q0, point_i);
    assert(~isempty(p.v_i), 'Error when looking for node on body %d', body_i_id)
end
p.idx = b.h_idx;

[b, isrigid] = h_get_body(sys, body_j_id);
if isrigid
    p.v_j = u_global_to_local(b, q0, point_j);
else
    p.v_j = u_get_node_id(b, q0, point_j);
    assert(~isempty(p.v_j), 'Error when looking for node on body %d', body_j_id)
end
p.idx = [p.idx, b.h_idx];

sys.force.tsd = [sys.force.tsd, p];

end

