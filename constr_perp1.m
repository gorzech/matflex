function [c, c_q, c_p, g] = constr_perp1(sys, d1_constr, q, h)
%CONSTR_PERP1 Constraint for perpendicular vectors

% This is a single constraint in the form vi'*vj
[b_i, isrigid] = h_get_body(sys, d1_constr.body_i_id);
if isrigid
    [v_i, cq_i, v_i_p, v_i_b_g] = d1_rigid(q, h, d1_constr.v_i, b_i);
else
    [v_i, cq_i, v_i_p, v_i_b_g] = d1_flex(q, h, d1_constr.v_i, d1_constr.psi_i, b_i);
end

[b_j, isrigid] = h_get_body(sys, d1_constr.body_j_id);
if isrigid
    [v_j, cq_j, v_j_p, v_j_b_g] = d1_rigid(q, h, d1_constr.v_j, b_j);
else
    [v_j, cq_j, v_j_p, v_j_b_g] = d1_flex(q, h, d1_constr.v_j, d1_constr.psi_j, b_j);
end

c = dot(v_i, v_j);
c_p = dot(v_i, v_j_p) + dot(v_j, v_i_p);
g = -dot(v_i, v_j_b_g) - dot(v_j, v_i_b_g) - 2 * dot(v_i_p, v_j_p);

c_q = zeros(1, sys.nh);
c_q(1, b_i.h_idx(4 : end)) = v_j' * cq_i;
c_q(1, b_j.h_idx(4 : end)) = v_i' * cq_j;
end

function [vi, cq_om, vi_p, vi_b_g] = d1_rigid(q, h, vi_local, body)
p = q(body.q_idx(4 : 7));
om = h(body.h_idx(4 : 6));
Ai = Rot(p);

vi = Ai * vi_local;
cq_om = -Ai * sksym(vi_local);
vi_p = cq_om * om;
vi_b_g = Ai * sksym(om) ^ 2 * vi_local;
end

function [vi, cq_om_qf, vi_p, vi_b_g] = d1_flex(q, h, vi_local, psi, body)
p = q(body.q_idx(4 : 7));
om = h(body.h_idx(4 : 6));
qf = q(body.q_idx(8 : end));
qf_p = h(body.h_idx(7 : end));
Ai = Rot(p);
% Transformation to extract rotation at attachment point
theta = psi * qf;
theta_p = psi * qf_p;
% Rotation matrix at attachment point - assumes small deformations
Af = eye(3) + sksym(theta);
Af_p = sksym(theta_p);

vi_def_local = Af * vi_local;
vi = Ai * vi_def_local;
cq_om_qf = [-Ai * sksym(vi_def_local), -Ai * sksym(vi_local) * psi];
vi_p = cq_om_qf * h(body.h_idx(4 : end));
vi_b_g = Ai * sksym(om) ^ 2 * vi_def_local + 2 * Ai * sksym(om) * (Af_p * vi_local);
end