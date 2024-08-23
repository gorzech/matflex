function Q = h_force_to_body(sys, q, F, body_id, v_i)
%H_FORCE_TO_BODY Helper to apply vector F to body b at v_i

[b, isrigid] = h_get_body(sys, body_id);
if isrigid
    L = c_dq_point_rigid(q, v_i, b);
else
    L = c_dq_point_flex(q, v_i, b);
end

Q = L' * F;

end

% Taken from constr_point
function L = c_dq_point_rigid(q, s_i, body)
    pp_i = q(body.q_idx(4 : 7));
    A_i = Rot(pp_i);
    L = [eye(3), -A_i * sksym(s_i)];
end

function L = c_dq_point_flex(q, location, body)
    % qn_idx - three coordinates for the qn vector - node location
    pp_i = q(body.q_idx(4 : 7));
    A_i = Rot(pp_i);
    s_i = location.m0 + location.m1 * q(body.q_idx(8 : end));
    L = [eye(3), -A_i * sksym(s_i), A_i * location.m1];
end
