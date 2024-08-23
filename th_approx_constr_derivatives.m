function [c_h, c_p, g] = th_approx_constr_derivatives(sys, cnstr_fun, q, h, t_base)
%TH_APPROX_CONSTR_DERIVATIVES Helper test function to approximate contraint
%derivatives
if nargin < 5
    t_base = 0;
end
    c_q = approxJacobian(@constraint_function, q);
    
    % c_p = c_q * h + c_t
    consrt_t_fun = @(t) cnstr_fun(sys, t, q, h);
    c_t = approxJacobian(consrt_t_fun, t_base);
    c_p = gamma_function(q) + c_t;

    cq_qp_dq = approxJacobian(@gamma_function, q);
    
    c_h = zeros(size(c_q, 1), length(h));
    ch_h_dh = c_h;
    % transform _dq terms to correspond _dh
    for ii = 1 : sys.nbodies
        b = sys.bodies(ii);
        c_h = transf_q_h(b, q, c_q, c_h);
        ch_h_dh = transf_q_h(b, q, cq_qp_dq, ch_h_dh);
    end
    for ii = 1 : sys.nfbodies
        b = sys.fbodies(ii);
        c_h = transf_q_h(b, q, c_q, c_h);
        ch_h_dh = transf_q_h(b, q, cq_qp_dq, ch_h_dh);
    end
    
    % Approximate c_tt using second order formula (after wiki)
    step = sqrt(eps) * 1000;
    c_tt = (consrt_t_fun(t_base + step) - 2 .* consrt_t_fun(t_base) ...
        + consrt_t_fun(t_base - step)) ./ step ^ 2;
    
    g = -ch_h_dh * h - c_tt;

    function c = constraint_function(qc)
        % Just in case
        qc = normalize_Euler_params(sys, qc);
        [c, ~, ~, ~] = cnstr_fun(sys, t_base, qc, h);
    end

    function g = gamma_function(qc) 
        qc = normalize_Euler_params(sys, qc);
        [~, b_cq, ~, ~] = cnstr_fun(sys, t_base, qc, h);
        g = b_cq * h; % Cq * qp
    end
end

function qc = normalize_Euler_params(sys, qc)
    for b = sys.bodies
        p_idx = b.q_idx(4 : 7);
        qc(p_idx) = qc(p_idx) ./ norm(qc(p_idx));
    end
    for b = sys.fbodies
        p_idx = b.q_idx(4 : 7);
        qc(p_idx) = qc(p_idx) ./ norm(qc(p_idx));
    end
end

function A_h = transf_q_h(b, q, A_q, A_h)
    A_h(:, b.h_idx(1 : 3)) = A_q(:, b.q_idx(1 : 3));
    LiT = Gep(q(b.q_idx(4 : 7)))';
    A_h(:, b.h_idx(4 : 6)) = 0.5 * A_q(:, b.q_idx(4 : 7)) * LiT;
    A_h(:, b.h_idx(7 : end)) = A_q(:, b.q_idx(8 : end));
end