function sys = h_sys_force_point(sys, force, location, q0, body_id)
%H_SYS_FORCE_POINT Function to generate base for point force

p.body_id = body_id;

[b, isrigid] = h_get_body(sys, body_id);
if isrigid
    p.v_i = u_global_to_local(b, q0, location);
else
    p.v_i = u_get_node_id(b, q0, location);
    assert(~isempty(p.v_i), 'Error when looking for node on body %d', body_id)
end

try
    F = force(0);
    assert(all(size(F) == [3, 1]), 'Incorrect force vector size')
catch
    error('Cannot evaluate force at time 0')
end

p.force = force;
p.idx = b.h_idx;

sys.force.point = [sys.force.point, p];

end

