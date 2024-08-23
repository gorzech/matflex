function yp = sys_ode(sys, t, y)
%SYS_ODE Generate derivative of the y vector for a multibody system
%   Code here uses hybrid approach, that is y = [q; h] where q = [r; p] and
%   h = [r_dot, omega]; where p are Euler parameters
%   According to Nikravesh, 13.2.2
   
    yp = zeros(sys.ny, 1);
    for b = sys.bodies
        reassign_coords(b)
    end
    for f = sys.fbodies
        reassign_coords(f)
    end
    M = mbd_mass(sys, y(1 : sys.nq));
    Q = mbd_force(sys, t, y(1 : sys.nq), y(sys.nq + 1 : end));
%     if sys.nconstr > 0
%         [constraints, c_dq, c_dt, gamma_v] = mbd_constr(sys, y(1 : sys.nq), y(sys.nq + 1 : end));
%         % compute inv(M) * Cq'
%         helper = M \ [(c_dq'), Q];
%         m_inv_c_dq_t = helper(:, 1 : end - 1);
%         
%         % compute Cq * inv(M) * Cq'
%         cq_m_cq = c_dq * m_inv_c_dq_t;
%         auc = helper(:, end); %M \ Q;
% 
%         % compute vector Cq*inv(M)*F - gamma + 2*alfa*Cp + beta*beta*C
%         lambda = c_dq * auc - gamma_v + sys.baumg_p1 * c_dt + sys.baumg_p2 * constraints;
%         % solve cq_m_cq * x = lambda to get x as true lambda
%         yp(sys.nq + 1 : end) = auc - m_inv_c_dq_t * (cq_m_cq \ lambda);
%     else
    if sys.nconstr > 0
        [constraints, c_dq, c_dt, gamma_v] = mbd_constr(sys, t, y(1 : sys.nq), y(sys.nq + 1 : end));
        M_cq = [M, c_dq'; c_dq, zeros(sys.nconstr)];
        rhs = [Q; gamma_v - sys.baumg_p1 * c_dt - sys.baumg_p2 * constraints];
        auc = M_cq \ rhs;

        yp(sys.nq + 1 : end) = auc(1 : sys.nh);
    else
        yp(sys.nq + 1 : end) = M \ Q;
    end
    function reassign_coords(body)
        qi = body.q_idx;
        hi = body.h_idx;
%         r = y(7 * ii - 6 : 7 * ii - 4);
        p = y(qi(4:7));
        p = p ./ norm(p);
        y(qi(4:7)) = p;
        rp = y(sys.nq + hi(1:3));
        om = y(sys.nq + hi(4:6));
        yp(qi(1:3)) = rp;
        % pp = ?
        L = Gep(p);
        yp(qi(4:7)) = 0.5 .* (L' * om);
        yp(qi(8:end)) = y(sys.nq + hi(7:end));
    end
end