function [q, elem_idx, node_idx] = straight_beam3d(element_length, n_elem)
%STRAIGHT_BEAM3D Generate coordinates and coordinate indices for straight beam
%   beam_length - total beam length 
%   n_elem - number of elements required
%   Output:
%   q - beam coordinates in initial state
%   elem_idx - q(elem_idx(:, i)) - gives coords of i-th element
    n_dof = 6;
    q = [linspace(0, element_length * n_elem, n_elem + 1)
        zeros(n_dof - 1, n_elem + 1)];
    q = q(:); % convert to column vector
    elem_idx = zeros(2 * n_dof, n_elem);
    row_elem_idx = 0 : n_dof : n_dof * n_elem - 1;
    for ii = 1 : 12
        elem_idx(ii, :) = ii + row_elem_idx;
    end
    
    node_idx = zeros(n_dof, n_elem + 1);
    row_node_idx = 0 : n_dof : n_dof * n_elem;
    for ii = 1 : 6
        node_idx(ii, :) = ii + row_node_idx;
    end
end

