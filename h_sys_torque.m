function sys = h_sys_torque(sys, torque, location, q0, body_id)
%H_SYS_TORQUE Function to generate base for torque applied to the body

p.body_id = body_id;
p.psi = [];

[b, isrigid] = h_get_body(sys, body_id);
if ~isrigid
    [~, p.psi] = u_get_node_id(b, q0, location);
    assert(~isempty(p.psi), 'Error when looking for node orientation on body %d', body_id)
end

try
    L = torque(0);
    assert(isequal(size(L), [3, 1]), 'Incorrect torque vector size')
catch
    error('Cannot evaluate force at time 0')
end

p.torque = torque;
p.idx = b.h_idx;

sys.force.torque = [sys.force.torque, p];

end

