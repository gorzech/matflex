function M = ffr_mass_sid(fbody, p, qf)
%FFR_MASS_SID Mass matrix based on the SID_FEM description

A = Rot(p);

m_tt = fbody.mass .* eye(3);

mCM = u_tay_mult(fbody.mCM, qf);
m_rt = sksym(A * mCM);

m_ft = fbody.Ct * A';

J = u_tay_mult(fbody.mmi, qf);
m_rr = [J(1), J(4), J(5)
    J(4), J(2), J(6)
    J(5), J(6), J(3)];

m_fr = u_tay_mult(fbody.Cr, qf);

m_ff = fbody.Me;

M = [m_tt, m_rt', m_ft'
    m_rt, m_rr, m_fr'
    m_ft, m_fr, m_ff];
end

