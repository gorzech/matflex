function Qv = ffr_quadratic_inertia_preintegrated(fbody, p, qf, om, qf_p)
%FFR_QUADRATIC_INERTIA_PREINTEGRATED Compute velocity dependent inertia
%forces for floating frame of reference formulation. Version based on
%minimum amount of matrix transformations and preintegrated FEM mass
%matrices
%   fbody - structure with information about flexible body
A = Rot(p);
oms = sksym(om);

Qv = zeros(fbody.nh, 1);

% Equations taken from the paper by Orzechowski et al. (2017)
Ij = fbody.m2 + fbody.m3 * qf;
Qv(1 : 3) = -A * oms * (oms * Ij + 2 * fbody.m3 * qf_p);

% Rotational part I7 & I8
m8 = fbody.m8;
m_rr = fbody.m7 + [m8(:, :, 1) * qf, m8(:, :, 2) * qf, m8(:, :, 3) * qf];
m_rr_p = [m8(:, :, 1) * qf_p, m8(:, :, 2) * qf_p, m8(:, :, 3) * qf_p];
% if use I9
if ~isempty(fbody.m9)
    m9 = fbody.m9;
    for ii = 1 : 3
        m_rr(ii, ii) = m_rr(ii, ii) + qf' * m9{ii} * qf;
        m_rr_p(ii, ii) = m_rr_p(ii, ii) + 2 .* qf' * m9{ii} * qf_p;
    end
    for ii = [4 : 6; 1, 1, 2; 2, 3, 3]
        m_ii = qf' * m9{ii(1)} * qf;
        m_ii_p = qf' * (m9{ii(1)} + m9{ii(1)}') * qf_p;
        m_rr(ii(2), ii(3)) = m_rr(ii(2), ii(3)) + m_ii;
        m_rr(ii(3), ii(2)) = m_rr(ii(3), ii(2)) + m_ii;
        m_rr_p(ii(2), ii(3)) = m_rr_p(ii(2), ii(3)) + m_ii_p;
        m_rr_p(ii(3), ii(2)) = m_rr_p(ii(3), ii(2)) + m_ii_p;
    end
end

% Rotational flexible part
m5T = fbody.m5T;
m_rf = fbody.m4 + ...
    [m5T(:, :, 1) * qf, m5T(:, :, 2) * qf, m5T(:, :, 3) * qf]';

Qv(4 : 6) = - oms * m_rr * om - m_rr_p * om - oms * m_rf * qf_p;

m_rf_om_d_om = [m5T(:, :, 1) * qf_p, m5T(:, :, 2) * qf_p, m5T(:, :, 3) * qf_p];
Qv(7 : end) = -2 .* m_rf_om_d_om * om;
for ii = 1 : 3
    Qv(7 : end) = Qv(7 : end) + 0.5 .* (om(ii) .* om' * m8(:, :, ii))';
end

if ~isempty(fbody.m9)
    m9_4_qf = (m9{4} + m9{4}') * qf;
    m9_5_qf = (m9{5} + m9{5}') * qf;
    m9_6_qf = (m9{6} + m9{6}') * qf;
    m9_col_1 = [ 2 * m9{1} * qf ...
        m9_4_qf ...
        m9_5_qf] * om;
    Qv(7 : end) = Qv(7 : end) + 0.5 .* om(1) * m9_col_1;
    m9_col_2 = [m9_4_qf ...
        2 * m9{2} * qf ...
        m9_6_qf] * om;
    Qv(7 : end) = Qv(7 : end) + 0.5 .* om(2) * m9_col_2;
    m9_col_3 = [m9_5_qf ...
        m9_6_qf ...
        2 * m9{3} * qf] * om;
    Qv(7 : end) = Qv(7 : end) + 0.5 .* om(3) * m9_col_3;
end

end

