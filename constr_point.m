function [c, c_q, c_p, g] = constr_point(sys, point_constr, q, h)
%CONSTR_POINT Constraint for common point

% Write for single constr
c_q = zeros(3, sys.nh);
[b, isrigid] = h_get_body(sys, point_constr.body_i_id);
if isrigid
    [c_i, c_q_i, c_p_i, g_i] = constr_point_rigid(q, h, point_constr.v_i, b);
else
    [c_i, c_q_i, c_p_i, g_i] = constr_point_flex(q, h, point_constr.v_i, b);
end
c_q(:, b.h_idx) = c_q_i;

[b, isrigid] = h_get_body(sys, point_constr.body_j_id);
if isrigid
    [c_j, c_q_j, c_p_j, g_j] = constr_point_rigid(q, h, point_constr.v_j, b);
else
    [c_j, c_q_j, c_p_j, g_j] = constr_point_flex(q, h, point_constr.v_j, b);
end
c_q(:, b.h_idx) = -c_q_j;
c = c_i - c_j;
c_p = c_p_i - c_p_j;
g = g_i - g_j;


end

function [c, c_q, c_p, g] = constr_point_rigid(q, h, s_i, body)
    R_i = q(body.q_idx(1 : 3));
    p_i = q(body.q_idx(4 : 7));
    om_i = h(body.h_idx(4 : 6));
    A_i = Rot(p_i);
    c = R_i + A_i * s_i;
    c_q = [eye(3), -A_i * sksym(s_i)];
    c_p = c_q * h(body.h_idx);
    g = -A_i * sksym(om_i) ^ 2 * s_i;
end

function [c, c_q, c_p, g] = constr_point_flex(q, h, location, body)
    R_i = q(body.q_idx(1 : 3));
    p_i = q(body.q_idx(4 : 7));
    om_i = h(body.h_idx(4 : 6));
    A_i = Rot(p_i);
    B2_i = location.m1;
    s_i = location.m0 + B2_i * q(body.q_idx(8 : end));
    c = R_i + A_i * s_i;
    c_q = [eye(3), -A_i * sksym(s_i), A_i * B2_i];
    c_p = c_q * h(body.h_idx);
    om_s = sksym(om_i);
    g = -A_i * om_s ^ 2 * s_i - 2 * A_i * om_s * B2_i * h(body.h_idx(7 : end));
end