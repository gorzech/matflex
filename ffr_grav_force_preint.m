function Qg = ffr_grav_force_preint(fbody, g, p, qf)
%FFR_GRAV_FORCE_PREINT Gravity force for the flexible body - preintegrated
%version    
    A = Rot(p);
    
    g_T_A = g' * A;
    m_tr = sksym(fbody.m2 + fbody.m3 * qf);
    rot_Fg = -(g_T_A * m_tr)';
    q_Fg = (g_T_A * fbody.m3)';
    Qg = [fbody.mass * g
        rot_Fg
        q_Fg];
end

