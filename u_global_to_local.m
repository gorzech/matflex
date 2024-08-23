function s = u_global_to_local(b, q0, s_global)
%U_GLOBAL_TO_LOCAL Transforms s_global to body b local coordinates
r = q0(b.q_idx(1 : 3));
p = q0(b.q_idx(4 : 7));
A = Rot(p);
s = A' * (s_global - r);
end

