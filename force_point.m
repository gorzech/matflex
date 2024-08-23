function Q = force_point(sys, point_force, t, q)
%FORCE_POINT Return value of the force point acting on the body

F = point_force.force(t);
Q = h_force_to_body(sys, q, F, point_force.body_id, point_force.v_i);
end