function M = beam3d_mass(rho, A, l, ~, Iy, Iz, lumped)
% Compute mass matrix M for 3D beam element based on 
% Przemieniecki, 1968, p. 294
% Similar to Beam4, but Ansys Beam4 has more features.
if nargin < 7
    lumped = false;
end

m = rho * A * l;
Jx = Iy + Iz; % According to ANSYS docs this should be a polar moment of inertia

if lumped
    % After ANSYS documentation
    M = diag([m/2, m/2, m/2, 0, 0, 0, m/2, m/2, m/2, 0, 0, 0]);
    return
end

m_diag = [1 / 3, 13 / 35 + 6 * Iz / (5 * A * l ^ 2), ...
    13 / 35 + 6 * Iy / (5 * A * l ^ 2), Jx / 3 / A, ...
    l ^ 2 / 105 + 2 * Iy / 15 / A, l ^ 2 / 105 + 2 * Iz / 15 / A];

M = diag([m_diag, m_diag]);
m_z_ry = -11 * l / 210 - Iy / (10 * A * l);
    
m_y_rz = 11 * l / 210 + Iz / (10 * A * l);

% coupling terms
m_yi_yj = 9 / 70 - 6 * Iz / (5 * A * l ^ 2);
m_zi_zj = 9 / 70 - 6 * Iy / (5 * A * l ^ 2);
m_rxi_rxj = Jx / 6 / A;
m_ryi_ryj = -l ^ 2 / 140 - Iy / 30 / A;
m_rzi_rzj = -l ^ 2 / 140 - Iz / 30 / A;

m_yj_rzi = 13 * l / 420 - Iz / (10 * A * l);
m_zj_ryi = -13 * l / 420 + Iy / (10 * A * l);

M(5, 3) = m_z_ry;
M(3, 5) = m_z_ry;
M(11, 9) = -m_z_ry;
M(9, 11) = -m_z_ry;

M(6, 2) = m_y_rz;
M(2, 6) = m_y_rz;
M(12, 8) = -m_y_rz;
M(8, 12) = -m_y_rz;

M(7, 1) = 1 / 6;
M(1, 7) = 1 / 6;

M(8, 2) = m_yi_yj;
M(2, 8) = m_yi_yj;
M(9, 3) = m_zi_zj;
M(3, 9) = m_zi_zj;
M(10, 4) = m_rxi_rxj;
M(4, 10) = m_rxi_rxj;
M(11, 5) = m_ryi_ryj;
M(5, 11) = m_ryi_ryj;
M(12, 6) = m_rzi_rzj;
M(6, 12) = m_rzi_rzj;

M(8, 6) = m_yj_rzi;
M(6, 8) = m_yj_rzi;
M(12, 2) = -m_yj_rzi;
M(2, 12) = -m_yj_rzi;
M(9, 5) = m_zj_ryi;
M(5, 9) = m_zj_ryi;
M(11, 3) = -m_zj_ryi;
M(3, 11) = -m_zj_ryi;

M = M .* m;