function Qg = ffr_grav_force(fbody, g, p, qf)
%FFR_GRAV_FORCE Gravity force for the flexible body
    elem_idx = fbody.elem_idx;
    
    S_int = beam3d_shape_fun_integrate(fbody.mat.rho, fbody.sec.A, fbody.elem_len);
    A = Rot(p);
    qn = fbody.q0 + fbody.B2 * qf;
    rot_Fg = zeros(3, 1);
    q_Fg = zeros(size(fbody.q0));
    
    
    g_T_A = g' * A;
    for ii = 1 : fbody.n_elem
        u = S_int * qn(elem_idx(:, ii));
        rot_Fg = rot_Fg - (g_T_A * sksym(u))';
        q_Fg(elem_idx(:, ii)) = q_Fg(elem_idx(:, ii)) + (g_T_A * S_int)';
    end
    Qg = [fbody.mass * g
        rot_Fg
        fbody.B2' * q_Fg];
end

