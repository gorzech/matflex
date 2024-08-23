function [c, c_q, c_p, g] = mbd_constr(sys, t, q, h)
%MBD_CONSTR Provide constraint vector and it derivatives for the system

c = zeros(sys.nconstr, 1);
c_p = c;
g = c;
c_q = zeros(sys.nconstr, sys.nh);

c_idx = 0;
for s = sys.constr.simple
    [c_i, c_q_i, c_p_i, g_i] = constr_simple(sys, s, q, h);
    set_constr_i;
end
for p = sys.constr.point
    [c_i, c_q_i, c_p_i, g_i] = constr_point(sys, p, q, h);
    set_constr_i;
end
for p = sys.constr.perp1
    [c_i, c_q_i, c_p_i, g_i] = constr_perp1(sys, p, q, h);
    set_constr_i;
end
for p = sys.constr.drive_point
    [c_i, c_q_i, c_p_i, g_i] = constr_drive_point(sys, p, t, q, h);
    set_constr_i;
end

    function set_constr_i()
        nc = length(c_i);
        c_i_idx = c_idx + (1 : nc);
        c(c_i_idx) = c_i;
        c_p(c_i_idx) = c_p_i;
        g(c_i_idx) = g_i;
        c_q(c_i_idx, :) = c_q_i;
        c_idx = c_idx + nc;
    end

end

