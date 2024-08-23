function b_pos = post_body_position(sys, T, Y, body_id, s, is_s_local)
%POST_BODY_POSITION Get position of the body point based on solution

if nargin < 5
    s = [0; 0; 0];
end 
if nargin < 6
    is_s_local = true;
end

[b, isrigid] = h_get_body(sys, body_id);

b_pos = zeros(3, length(T));
if isrigid
    if ~is_s_local
        [r, A] = get_r_A(b, Y, 1);
        s_loc = A' * (s - r);
    else
        s_loc = s;
    end
    for ii = 1 : length(T)
        [r, A] = get_r_A(b, Y, ii);
        b_pos(:, ii) = r + A * s_loc;
    end
else
    if is_s_local
        [r, A] = get_r_A(b, Y, 1);
        s_glob = r + A * s;
    else
        s_glob = s;
    end
    v = u_get_node_id(b, Y(1, 1 : sys.nq)', s_glob);
    for ii = 1 : length(T)
        [r, A] = get_r_A(b, Y, ii);
        qf = Y(ii, b.q_idx(8 : end))';
        qn = v.m0 + v.m1 * qf;
        b_pos(:, ii) = r + A * qn;
    end
end

end

function [r, A] = get_r_A(body, Y, idx)
    r = Y(idx, body.q_idx(1 : 3))';
    p = Y(idx, body.q_idx(4 : 7))';
    A = Rot(p);
end