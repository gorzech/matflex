function Qg = ffr_grav_force_sid(fbody, g, p, qf)
%FFR_GRAV_FORCE_PREINT Gravity force for the flexible body - preintegrated
%version    
    A = Rot(p);
    
    A_T_g = A' * g;
    mCM = u_tay_mult(fbody.mCM, qf);
    m_rt = sksym(mCM);
    
    rot_Fg = m_rt * A_T_g;
    q_Fg = fbody.Ct * A_T_g;
    Qg = [fbody.mass * g
        rot_Fg
        q_Fg];
end

