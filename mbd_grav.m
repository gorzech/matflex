function g = mbd_grav(sys, q)
%MBD_GRAV Computes gravity vector
    % Assumes sys contains nh, nbodies and grav
    g = zeros(sys.nh, 1);
    for b = sys.bodies
        g(b.h_idx(1:3)) = b.mass .* sys.grav;
    end
    for f = sys.fbodies
        g(f.h_idx) = f.ffr_grav(f, sys.grav, q(f.q_idx(4:7)), q(f.q_idx(8:end)));
    end
end

