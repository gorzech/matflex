function qv = mbd_quadratic_inertia(sys, q, h)
%MBD_QUADRATIC_INERTIA Returns so-called quadratic velocity vector
    
    qv = zeros(sys.nh, 1);
    for b = sys.bodies
        om = h(b.h_idx(4 : 6));
        qv(b.h_idx(4 : 6)) = cross3(om, b.Ic * om);
    end
    for f = sys.fbodies
        qv(f.h_idx) = ...
            f.ffr_quad_inertia(f, q(f.q_idx(4:7)), q(f.q_idx(8:end)), ...
            h(f.h_idx(4 : 6)), h(f.h_idx(7:end)));
    end
end
