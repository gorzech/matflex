function Q = mbd_force(sys, t, q, h)
%MBD_FORCE Get total system forces
    if sys.usequadratic
        Q = -mbd_quadratic_inertia(sys, q, h);
    else
        Q = zeros(size(h));
    end
    if sys.isgrav
        Q = Q + mbd_grav(sys, q);
    end
    % add elastic forces (if any)
    for f = sys.fbodies
        if isfield(f, 'De')
            Q(f.h_idx(7:end)) = Q(f.h_idx(7:end)) - f.Kff * q(f.q_idx(8:end)) ... 
                - f.De * h(f.h_idx(7:end));
        else
            Q(f.h_idx(7:end)) = Q(f.h_idx(7:end)) - f.Kff * q(f.q_idx(8:end));
        end
    end
    for p = sys.force.point
        Q(p.idx) = Q(p.idx) + force_point(sys, p, t, q);
    end
    for p = sys.force.tsd
        Q(p.idx) = Q(p.idx) + force_tsd(sys, p, q, h);
    end
    for p = sys.force.torque
        Q(p.idx) = Q(p.idx) + force_torque(sys, p, t, q);
    end
end

