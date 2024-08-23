function [c, c_q, c_p, g] = constr_drive_point(sys, point_drv_constr, t, q, h)
%CONSTR_POINT Constraint for common point

% Write for single constr
c_q = zeros(point_drv_constr.nc, sys.nh);
[b, isrigid] = h_get_body(sys, point_drv_constr.body_id);
if isrigid
    [c_i, c_q_i, c_p_i, g_i] = constr_point_rigid(q, h, point_drv_constr.v, b);
else
    [c_i, c_q_i, c_p_i, g_i] = constr_point_flex(q, h, point_drv_constr.v, b);
end
c_q(:, b.h_idx) = c_q_i(point_drv_constr.idx, :);

c = c_i(point_drv_constr.idx) - point_drv_constr.constr_fun(t);
c_p = c_p_i(point_drv_constr.idx) - point_drv_constr.constr_fun_dt(t);
g = g_i(point_drv_constr.idx) + point_drv_constr.constr_fun_d2t(t);

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