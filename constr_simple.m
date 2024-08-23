function [c, c_q, c_p, g] = constr_simple(sys, simple_constr, q, h)
%CNSTR_SIMPLE Simple constraint - removes certain DOFs globally
%   Return constraint vector c, Jacobian c_q, derivative c_p and gamma
%   vector g

s = simple_constr;
nc = s.nc;

c = q(s.q_idx) - s.q0;
g = zeros(nc, 1);

c_q = zeros(nc, sys.nh);
r_idx = (1 : s.nc - 3);
p_idx = (s.nc - 2 : s.nc);
c_q(r_idx, s.h_idx(r_idx)) = eye(s.nc - 3);
c_q(p_idx, s.h_idx(p_idx)) = s.G0_2;
c_p = h(s.h_idx);
c_p(p_idx) = s.G0_2 * c_p(p_idx);
% g is equal to zero

