function [location, orientation] = u_get_node_id(b, q0, s_global)
%U_GET_NODE_ID how to get to that point of flexible body using 
% location = m0 + m1 * qf
% orientation * qf - rotational DOFs

s = u_global_to_local(b, q0, s_global);
% Need to find s in q0
location = [];
orientation = [];
if isfield(b, 'q0') % not a SID fbody
    for ii = 1 : 6 : length(b.q0)
        l = norm(s - b.q0(ii : ii + 2));
        if l < 1e-5
            location = struct('m0', b.q0(ii : ii + 2), 'm1', b.B2(ii : ii + 2, :));
            orientation = b.B2(ii + 3 : ii + 5, :);
            break;
        end
    end
else % SID body
    for ii = 1 : b.n_node
        l = norm(s - b.node(ii).origin.m0);
        if l < 1e-5
            location = b.node(ii).origin;
            orientation = b.node(ii).psi;
            break;
        end
    end
end
end

