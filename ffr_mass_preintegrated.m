function M = ffr_mass_preintegrated(fbody, p, qf)
%FFR_MASS_PREINTEGRATED Procedure to get mass matrix based on direct
%use of the FEM mass matrix 
%   fbody contains mass integrals
%   
%   p - orientation vector - Eueler parameters
%   qf - flexible parameters. Originals as for fem analysis

A = Rot(p); % Rotation matrix

% Translational part
m_tt = fbody.mass .* speye(3); 

% Rotational part I7 & I8
m8 = fbody.m8;
m_rr = fbody.m7 + [m8(:, :, 1) * qf, m8(:, :, 2) * qf, m8(:, :, 3) * qf];
% if use I9
if ~isempty(fbody.m9)
    m9 = fbody.m9;
    for ii = 1 : 3
        m_rr(ii, ii) = m_rr(ii, ii) + qf' * m9{ii} * qf;
    end
    for ii = [4 : 6; 1, 1, 2; 2, 3, 3]
        m_ii = qf' * m9{ii(1)} * qf; % Same for "transposed version"
        m_rr(ii(2), ii(3)) = m_rr(ii(2), ii(3)) + m_ii;
        m_rr(ii(3), ii(2)) = m_rr(ii(3), ii(2)) + m_ii;
    end
end

% Flexible component
m_ff = fbody.m_ff;

% Tranlational-rotational component
m_tr = sksym(fbody.m2 + fbody.m3 * qf);
m_tr = -A * m_tr;

% Translational flexible part
m_tf = A * fbody.m3;

% Rotational flexible part
m5T = fbody.m5T;
m_rf = fbody.m4 + ...
    [m5T(:, :, 1) * qf, m5T(:, :, 2) * qf, m5T(:, :, 3) * qf]';

M = [m_tt, m_tr, m_tf
    m_tr', m_rr, m_rf
    m_tf', m_rf', m_ff];
end