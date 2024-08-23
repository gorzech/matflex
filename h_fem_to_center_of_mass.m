function [q0, mass] = h_fem_to_center_of_mass(M, q0, node_idx)
%H_FEM_TO_CENTER_OF_MASS Set q0 w.r.t. flexible body center of mass
    n_node = size(node_idx, 2);

    % Find center of mass
    phi_rigid = zeros(6, length(q0));
    for ii = 1 : n_node
        phi_rigid(:, node_idx(:, ii)) = [eye(3), zeros(3)
            sksym(q0(node_idx(1 : 3, ii))), eye(3)];
    end
    M_rigid = phi_rigid * M * phi_rigid';
    mass = M_rigid(1, 1);
    % est = [0 -e(3) e(2)
    %   e(3) 0 -e(1)
    %   -e(2) e(1) 0];
    % Position of the center of mass from "inverse" sksym 
    R = [M_rigid(3, 5); M_rigid(1, 6); M_rigid(2, 4)] ./ mass;
    for ii = 1 : n_node
        % shift position
        q0(node_idx(1 : 3, ii)) = q0(node_idx(1 : 3, ii)) + R;
    end
end