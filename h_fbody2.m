function fbody = h_fbody2(q0, M, K, V, node_idx, move_to_cm)
%H_FBODY Heleper to generate flexible body - accepts finite element mass
%and stiffness matrix - also requires transformation V and node_idx

if nargin < 6
    move_to_cm = false;
end

n_node = size(node_idx, 2);

if move_to_cm
    q0 = h_fem_to_center_of_mass(M, q0, node_idx);
end
fbody.q0 = q0;
fbody.n_node = n_node;
fbody.node_idx = node_idx;

M_preint_mult = zeros(3, length(q0));
for ii = 1 : n_node
    M_preint_mult(:, node_idx(1 : 3, ii)) = eye(3);
end
% Effectively this is preintegrated S for whole body w.r.t. rho dV
S_int = M_preint_mult * M;
mass = 0;
for ii = 1 : n_node
    mass = mass + S_int(1, node_idx(1, ii));
end
fbody.mass = mass;

fbody.B2 = V;

fbody.m_ff = V' * M * V;

fbody.nflex = size(V, 2);
fbody.nh = 6 + fbody.nflex;
fbody.nq = fbody.nh + 1;

fbody.Kff = V' * K * V;

[fbody.m7, fbody.m8, fbody.m9] = m_rr_sum(node_idx, n_node, q0, S_int, V, mass);

fbody.m2 = S_int * q0;
fbody.m3 = S_int * V;

[fbody.m4, fbody.m5T] = m_rf_sum(node_idx, n_node, q0, S_int, fbody.m9, V);

end

% Mass integrals for the m_rf part
function [m4, m5T] = m_rf_sum(node_idx, n_node, q0, S_int, m9, B2)
    nq = size(m9{1}, 1);
    m4 = zeros(3, length(q0));
    for ii = 1 : n_node
        u0 = q0(node_idx(1:3, ii));
        u0s = sksym(u0);
        m4(:, node_idx(:, ii)) = u0s * S_int(:, node_idx(:, ii));
    end
    m4 = m4 * B2;
    % According to the paper this is based on I9 parts!
    m5T = zeros(nq, nq, 3); % transposition
    m5T(:, :, 1) = m9{6} - m9{6}';
    m5T(:, :, 2) = m9{5}' - m9{5};
    m5T(:, :, 3) = m9{4} - m9{4}';
end

% Mass integrals for the m_rr part
function [m7, m8, m9] = m_rr_sum(node_idx, n_node, q0, S_int, B2, mass)
    nq = size(B2, 2);
    m7 = zeros(3);
    m8 = zeros([3, nq, 3]);
    m9 = cell(1, 6); 
    m9(:) = {sparse(nq, nq)};
    S_int_sq2 = sign(S_int) .* sqrt(abs(S_int));
    % m_rr = sum of sksym(qn)' * sksym(qn)
    % and m8 - qn = q0 + B2 * qf;
    for ii = 1 : n_node
        node_mass = S_int(1, node_idx(1, ii));
        u0 = q0(node_idx(1:3, ii));
        us0 = sksym(u0);
        m7 = m7 + us0' * us0 .* node_mass;
        % Collect all q0 * B2 parts
        B2_i = S_int(:, node_idx(:, ii)) * B2(node_idx(:, ii), :);
        for jj = 1 : nq
            B2_ij_s = sksym(B2_i(:, jj));
            ts = us0' * B2_ij_s + B2_ij_s' * us0;
%             ts = ts .* node_mass;
            for kk = 1 : 3
                m8(:, jj, kk) = m8(:, jj, kk) + ts(:, kk);
            end
        end
        % now m9 - B2' * B2
        % Base on symbolic sksym(u)' * sksym(u) version
        % m9 contains entries 11, 22, 33, 12, 13, 23
        B2_i = S_int_sq2(:, node_idx(:, ii)) * B2(node_idx(:, ii), :);
        m9{1} = m9{1} + (B2_i(2, :)' * B2_i(2, :) + B2_i(3, :)' * B2_i(3, :));
        m9{2} = m9{2} + (B2_i(1, :)' * B2_i(1, :) + B2_i(3, :)' * B2_i(3, :));
        m9{3} = m9{3} + (B2_i(1, :)' * B2_i(1, :) + B2_i(2, :)' * B2_i(2, :));
        m9{4} = m9{4} - (B2_i(1, :)' * B2_i(2, :));
        m9{5} = m9{5} - (B2_i(1, :)' * B2_i(3, :));
        m9{6} = m9{6} - (B2_i(2, :)' * B2_i(3, :));
    end
    % Check for a beam - if any Iii is close to zero
    m7_diag = diag(m7);
    idx = find(m7_diag < 1e-6 * sum(m7_diag));
    if length(idx) > 1
        warning("There are large differences in flexible body inertia tensor - please check")
    elseif length(idx) == 1
        % seems we have beam here
        % Get beam length
        cmin = q0(node_idx(idx, 1));
        cmax = cmin;
        for ii = 2 : n_node
            ccurr = q0(node_idx(idx, ii));
            cmin = min(cmin, ccurr);
            cmax = max(cmax, ccurr);
        end
        l = cmax - cmin;
        others = [1 : idx - 1, idx + 1 : 3];
        % approximate missing I assuming cylinder
        m7(idx, idx) = m7(others(1), others(1)) + m7(others(2), others(2))...
            - mass * l ^ 2 / 6;
    end
end
