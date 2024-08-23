function [b, isrigid] = h_get_body(sys, body_id)
%H_GET_BODY Small helper to return body base on its id
if body_id <= sys.nbodies
    b = sys.bodies(body_id);
    isrigid = true;
else
    b = sys.fbodies(body_id - sys.nbodies);
    isrigid = false;
end

