function Q = force_tsd(sys, tsd_force, q, h)
%FORCE_TSD Return translational spring damper force between two bodies

[b, isrigid] = h_get_body(sys, tsd_force.body_i_id);
if isrigid
    [l_i, L_i, l_p_i] = constr_point_rigid(q, h, tsd_force.v_i, b);
else
    [l_i, L_i, l_p_i] = constr_point_flex(q, h, tsd_force.v_i, b);
end

[b, isrigid] = h_get_body(sys, tsd_force.body_j_id);
if isrigid
    [l_j, L_j, l_p_j] = constr_point_rigid(q, h, tsd_force.v_j, b);
else
    [l_j, L_j, l_p_j] = constr_point_flex(q, h, tsd_force.v_j, b);
end
l_vector = l_j - l_i;
dist = norm(l_vector);

l_p_vector = l_p_j - l_p_i;
dist_p = (l_vector' * l_p_vector) ./ dist;

f_value = tsd_force.k * (dist - tsd_force.l0) + tsd_force.d * dist_p;

F = (f_value ./ dist) * l_vector;

Q = [L_i' * F; -L_j' * F];

end

% Taken from constr_point
function [l, L, l_p] = constr_point_rigid(q, h, s_i, body)
    R_i = q(body.q_idx(1 : 3));
    pp_i = q(body.q_idx(4 : 7));
    A_i = Rot(pp_i);
    l = R_i + A_i * s_i;
    L = [eye(3), -A_i * sksym(s_i)];
    l_p = L * h(body.h_idx);
end

function [l, L, l_p] = constr_point_flex(q, h, location, body)
    % qn_idx - three coordinates for the qn vector - node location
    R_i = q(body.q_idx(1 : 3));
    pp_i = q(body.q_idx(4 : 7));
    A_i = Rot(pp_i);
    s_i = location.m0 + location.m1 * q(body.q_idx(8 : end));
    l = R_i + A_i * s_i;
    L = [eye(3), -A_i * sksym(s_i), A_i * location.m1];
    l_p = L * h(body.h_idx);
end
