function sys = h_sys_constr_drive_point(sys, location, q0, body_id, pos_idx, pos_fun, pos_fun_dt, pos_fun_d2t)
%H_SYS_CONSTR_DRIVE_POINT Helper function for adding drive constraint
%parameters for simple constraint that fixes point location as a function 
%of time

p.nc = length(pos_idx);
p.body_id = body_id;
p.idx = pos_idx;
p.constr_fun = pos_fun;
p.constr_fun_dt = pos_fun_dt;
p.constr_fun_d2t = pos_fun_d2t;

[b, isrigid] = h_get_body(sys, body_id);
if isrigid
    p.v = u_global_to_local(b, q0, location);
else
    p.v = u_get_node_id(b, q0, location);
    assert(~isempty(p.v), 'Error when looking for node on body %d', body_id)
end

sys.nconstr = sys.nconstr + p.nc;
sys.constr.drive_point = [sys.constr.drive_point, p];

end