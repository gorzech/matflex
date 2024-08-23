function Qv = ffr_quadratic_inertia_sid(fbody, p, qf, om, qf_p)
%FFR_QUADRATIC_INERTIA_PREINTEGRATED Compute velocity dependent inertia
%forces for floating frame of reference formulation. Version based on
%SID flexible body description
%   fbody - structure with information about flexible body
A = Rot(p);
oms = sksym(om);
nq = fbody.nflex;

Qv = zeros(fbody.nh, 1);

mCM = u_tay_mult(fbody.mCM, qf);
Qv(1 : 3) = -A * oms * (oms * mCM + 2 * fbody.Ct' * qf_p);

J = u_tay_mult(fbody.mmi, qf);
m_rr = [J(1), J(4), J(5)
    J(4), J(2), J(6)
    J(5), J(6), J(3)];

Gr = [fbody.Gr(:, 1 : nq) * qf_p, ...
    fbody.Gr(:, nq + 1 : 2 * nq) * qf_p, ...
    fbody.Gr(:, 2 * nq + 1 : 3 * nq) * qf_p] * om;

Qv(4 : 6) = - oms * (m_rr * om) - Gr;

Ge = [fbody.Ge(:, 1 : nq) * qf_p, ...
    fbody.Ge(:, nq + 1 : 2 * nq) * qf_p, ...
    fbody.Ge(:, 2 * nq + 1 : 3 * nq) * qf_p] * om;

Oe = u_tay_mult(fbody.Oe, qf);

OM = [om .^ 2; om(1) * om(2); om(2) * om(3); om(1) * om(3)];

Qv(7 : end) = -Ge - Oe * OM;

end

