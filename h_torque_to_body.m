function Q = h_torque_to_body(sys, q, L, body_id, psi)
%H_TORQUE_TO_BODY Helper to apply torque L to body b with the help of psi

[b, isrigid] = h_get_body(sys, body_id);
pp_i = q(b.q_idx(4 : 7));
A_i = Rot(pp_i)';
L_i = A_i * L;

if isrigid
    Q = [zeros(3,1); L_i];
else
    Q = [zeros(3,1); L_i; psi' * L_i];
end

end
