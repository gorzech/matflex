function fbody = h_fbody3_J(q0, M, K, V, node_idx)
%H_FBODY Heleper to generate flexible body - accepts finite element mass
%and stiffness matrix - also requires transformation V and node_idx
%Version that tries to mimic SID generation procedure - thus includes
%rotations

n_node = size(node_idx, 2);

[q0, mass] = h_fem_to_center_of_mass(M, q0, node_idx);

fbody.mass = mass;
fbody.q0 = q0;
fbody.n_node = n_node;
fbody.node_idx = node_idx;
fbody.B2 = V;

phi_rigid = zeros(6, length(q0));
for ii = 1 : n_node
    phi_rigid(:, node_idx(:, ii)) = [eye(3), zeros(3)
        sksym(q0(node_idx(1 : 3, ii))), eye(3)];
end
% Effectively this is preintegrated S for whole body w.r.t. rho dV
S_int = phi_rigid * M;
M_rigid = S_int * phi_rigid';

fbody.m_ff = V' * M * V;

fbody.nflex = size(V, 2);
fbody.nh = 6 + fbody.nflex;
fbody.nq = fbody.nh + 1;

fbody.Kff = V' * K * V;

fbody.m7 = M_rigid(4:6, 4:6);

% m8_v1 = m_rr_qf_1(node_idx, n_node, q0, S_int, V);
m8_v2 = m_rr_qf_2(node_idx, n_node, q0, S_int, V, phi_rigid); % faster

% fprintf('\n\t%g - norm between m8\n\n', norm(m8_v1(:) - m8_v2(:)))
fbody.m8 = m8_v2;

% fbody.m9 = m_rr_qf_qf(node_idx, n_node, S_int, V);
fbody.m9 = m_rr_qf_qf_2(node_idx, n_node, S_int, V);

% For m_tr
fbody.m2 = zeros(3, 1); %S_int(1:3, :) * q0;
fbody.m3 = S_int(1:3, :) * V; % Also (mainly?) m_tf

% Const part of m_rf
fbody.m4 = S_int(4:6, :) * V;
fbody.m5T = m_rf_sum(fbody.m9);

end

% Mass integrals for the m_rf part
function m5T = m_rf_sum(m9)
    nq = size(m9{1}, 1);
    % According to the paper this is based on I9 parts!
    m5T = zeros(nq, nq, 3); % transposition
    m5T(:, :, 1) = m9{6} - m9{6}';
    m5T(:, :, 2) = m9{5}' - m9{5};
    m5T(:, :, 3) = m9{4} - m9{4}';
end

% % Mass integrals for the m_rr part
% function m8 = m_rr_qf_1(node_idx, n_node, q0, S_int, V)
%     nq = size(V, 2);
%     m8 = zeros([3, nq, 3]);
%     % m_rr = sum of sksym(qn)' * sksym(qn)
%     % and m8 - qn = q0 + B2 * qf;
%     for ii = 1 : n_node
%         u0 = q0(node_idx(1:3, ii));
%         us0 = sksym(u0);
%         % Collect all q0 * B2 parts
%         B2_i = S_int(:, node_idx(:, ii)) * V(node_idx(:, ii), :);
%         for jj = 1 : nq
%             B2_ij_s = sksym(B2_i(:, jj));
%             ts = us0' * B2_ij_s + B2_ij_s' * us0;
% %             ts = ts .* node_mass;
%             for kk = 1 : 3
%                 m8(:, jj, kk) = m8(:, jj, kk) + ts(:, kk);
%             end
%         end
%     end
% end

function m8 = m_rr_qf_2(node_idx, n_node, q0, S_int, V, phi_rigid)
    nq = size(V, 2);
    qf_preint = zeros(3, length(q0), nq);
    for ii = 1 : n_node
        inode = node_idx(:, ii);
        qf_i = S_int(:, inode) * V(inode, :);
        for jj = 1 : nq
            qf_preint(:, inode, jj) = [...
                sksym(qf_i(1 : 3, jj)), diag( qf_i(4 : 6, jj) )]; %diag( B2_i(4 : 6, jj) )];
        end
    end
    m8 = zeros([3, nq, 3]);
    for jj = 1 : nq
        m8_jj = phi_rigid(4:6, :) * qf_preint(:, :, jj)';
        m8(:, jj, :) = m8_jj + m8_jj';
    end
end

% function m8_mmi = m_rr_qf_3(node_idx, n_node, q0, M, V)
%     qn0 = zeros(6, length(q0));
%     for ii = 1 : n_node
%         inode = node_idx(:, ii);
%     %     q0s = sksym(q0(inode(1 : 3)));
%     %     J_q0 = q0s' * q0s;
%         qn0(1, inode([2, 3])) = 2 .* q0(inode([2, 3]));
%         qn0(2, inode([1, 3])) = 2 .* q0(inode([1, 3]));
%         qn0(3, inode([1, 2])) = 2 .* q0(inode([1, 2]));
%         qn0(4, inode([1, 2])) = -q0(inode([2, 1]));
%         qn0(4, inode(6)) = -2;
%         qn0(5, inode([1, 3])) = -q0(inode([3, 1]));
%         qn0(5, inode(5)) = -2;
%         qn0(6, inode([2, 3])) = -q0(inode([3, 2]));
%         qn0(6, inode(4)) = -2;
%     end
%     % 
%     m8_mmi = qn0 * M * V;
%     m8_mmi(abs(m8_mmi) < 2e-11) = 0;
% end

% function m9 = m_rr_qf_qf(node_idx, n_node, S_int, V)
%     nq = size(V, 2);
%     m9 = cell(1, 6); 
%     m9(:) = {sparse(nq, nq)};
%     S_int_sq2 = sign(S_int) .* sqrt(abs(S_int));
%     % m_rr = sum of sksym(qn)' * sksym(qn)
%     % and m8 - qn = q0 + B2 * qf;
%     for ii = 1 : n_node
%         % now m9 - B2' * B2
%         % Base on symbolic sksym(u)' * sksym(u) version
%         % m9 contains entries 11, 22, 33, 12, 13, 23
%         B2_i = S_int_sq2(:, node_idx(:, ii)) * V(node_idx(:, ii), :);
%         m9{1} = m9{1} + (B2_i(2, :)' * B2_i(2, :) + B2_i(3, :)' * B2_i(3, :));
%         m9{2} = m9{2} + (B2_i(1, :)' * B2_i(1, :) + B2_i(3, :)' * B2_i(3, :));
%         m9{3} = m9{3} + (B2_i(1, :)' * B2_i(1, :) + B2_i(2, :)' * B2_i(2, :));
%         m9{4} = m9{4} - (B2_i(1, :)' * B2_i(2, :));
%         m9{5} = m9{5} - (B2_i(1, :)' * B2_i(3, :));
%         m9{6} = m9{6} - (B2_i(2, :)' * B2_i(3, :));
%     end
% end

function m9 = m_rr_qf_qf_2(node_idx, n_node, S_int, V)
    nq = size(V, 2);
    m9 = cell(1, 6); 
    m9(:) = {sparse(nq, nq)};
    % m_rr = sum of sksym(qn)' * sksym(qn)
    % and m8 - qn = q0 + B2 * qf;
    for jj = 1 : n_node
        inode = node_idx(1:3, jj);
        % now m9 - B2' * B2
        % Base on symbolic sksym(u)' * sksym(u) version
        % m9 contains entries 11, 22, 33, 12, 13, 23
        B2_i = S_int(1 : length(inode), inode) * V(inode, :);
        B2_j = B2_i ./ S_int(1, inode(1));
%             B2_j = phi_rigid(1 : length(inode), inode) * V(inode, :);
        m9{1} = m9{1} + (B2_i(2, :)' * B2_j(2, :) + B2_i(3, :)' * B2_j(3, :));
        m9{2} = m9{2} + (B2_i(1, :)' * B2_j(1, :) + B2_i(3, :)' * B2_j(3, :));
        m9{3} = m9{3} + (B2_i(1, :)' * B2_j(1, :) + B2_i(2, :)' * B2_j(2, :));
        m9{4} = m9{4} - (B2_i(1, :)' * B2_j(2, :));
        m9{5} = m9{5} - (B2_i(1, :)' * B2_j(3, :));
        m9{6} = m9{6} - (B2_i(2, :)' * B2_j(3, :));
    end
end
