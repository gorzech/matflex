function M = ffr_mass_integrate(fbody, p, qf)
%FFR_MASS_INTEGRATE Procedure to get mass matrix based on direct
%integration of the shape function (for most parts)
%   fbody contains (at lest)
%   q0 - initial coordinates
%   B2 - reference conditions
%   elem_idx - idices of the elements. Needed for integration
%   mat - material properties (rho)
%   sec - section props (Area, thky, thkz)
%   elem_len - length of the elements
%   n_elem - number of elements
%   fem_mass - FEM mass matrix
%   
%   p - orientation vector - Eueler parameters
%   qf - flexible parameters. Originals as for fem analysis

B2 = fbody.B2;
elem_idx = fbody.elem_idx;
le = fbody.elem_len;
n_elem = fbody.n_elem;

vint = @(f, rx, ry, rz) volumetric_integral(f, fbody.mat, fbody.sec, le, rx, ry, rz);
qn = fbody.q0 + B2 * qf;
A = Rot(p); % Rotation matrix

% Translational part
m_tt = fbody.mass .* speye(3); 
% Rotational part
int_mrr = @(xi, eta, zeta) m_rr_int(xi, eta, zeta, le, elem_idx, n_elem, qn);

m_rr = vint(int_mrr, 4, 2, 2);
if m_rr(1, 1) < sqrt(eps(m_rr(2, 2))) % beam?!
    m_rr(1, 1) = m_rr(2, 2) + m_rr(3, 3);
%     m_rr(1, 1) = sqrt(eps(m_rr(2, 2)));
end

% Flexible component
m_ff = fbody.m_ff;

% Tranlational-rotational component
int_mtr = @(xi, eta, zeta) m_tr_int(xi, eta, zeta, le, elem_idx, n_elem, qn);
m_tr = vint(int_mtr, 2, 1, 1);
m_tr = -A * m_tr;

% Translational flexible part
int_mtf = @(xi, eta, zeta) m_tf_int(xi, eta, zeta, le, elem_idx, n_elem, length(qn));
m_tf = vint(int_mtf, 2, 1, 1);
m_tf = A * m_tf * B2;

% Rotational flexible part
int_mrf = @(xi, eta, zeta) m_rf_int(xi, eta, zeta, le, elem_idx, n_elem, qn);
m_rf = vint(int_mrf, 4, 2, 2) * B2;

M = [m_tt, m_tr, m_tf
    m_tr', m_rr, m_rf
    m_tf', m_rf', m_ff];

end

% integrant of the m_rf part
function m_rf = m_rf_int(xi, eta, zeta, le, elem_idx, n_elem, qn)
    m_rf = zeros(3, length(qn));
    S = beam3d_shape_fun(xi, eta, zeta, le);
    for ii = 1 : n_elem
        u = u_ij(S, ii, elem_idx, qn);
        us = sksym(u);
        m_rf(:, elem_idx(:, ii)) = m_rf(:, elem_idx(:, ii)) + us * S;
    end
end

% integrant of the m_tf part
function m_tf = m_tf_int(xi, eta, zeta, le, elem_idx, n_elem, n_qn)
    m_tf = zeros(3, n_qn);
    S = beam3d_shape_fun(xi, eta, zeta, le);
    for ii = 1 : n_elem
        m_tf(:, elem_idx(:, ii)) = m_tf(:, elem_idx(:, ii)) + S;
    end
end

% integrant of the m_tr part
function m_tr = m_tr_int(xi, eta, zeta, le, elem_idx, n_elem, qn)
    m_tr = zeros(3);
    S = beam3d_shape_fun(xi, eta, zeta, le);
    for ii = 1 : n_elem
        u = u_ij(S, ii, elem_idx, qn);
        us = sksym(u);
        m_tr = m_tr + us;
    end
end

% integrant of the m_rr part
function m_rr = m_rr_int(xi, eta, zeta, le, elem_idx, n_elem, qn)
    m_rr = zeros(3);
    S = beam3d_shape_fun(xi, eta, zeta, le);
    for ii = 1 : n_elem
        u = u_ij(S, ii, elem_idx, qn);
        us = sksym(u);
        m_rr = m_rr + us' * us;
    end
end

function u = u_ij(S, i_elem, elem_idx, qn)
    u = S * qn(elem_idx(:, i_elem));
end