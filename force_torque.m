function Q = force_torque(sys, torque, t, q)
%FORCE_TORQUE Return value of the torque acting on the body

L = torque.torque(t);
Q = h_torque_to_body(sys, q, L, torque.body_id, torque.psi);
end