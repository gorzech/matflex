function M = mbd_mass(sys, q)
%MBD_MASS Computes mass matrix of the multibody system
    if sys.nfbodies > 0
        Mf = cell(1, sys.nfbodies);
        for ii = 1 : sys.nfbodies
            f = sys.fbodies(ii);
            Mf{ii} = f.ffr_mass(f, q(f.q_idx(4:7)), q(f.q_idx(8:end)));
        end
        M = blkdiag(sys.m_rb, Mf{:});
    else
        M = sys.m_rb;
    end
end

