function Qv = ffr_quadratic_inertia_simple_integrate(fbody, p, qf, om, qf_p)
%FFR_QUADRATIC_INERTIA_SIMPLE_INTEGRATE Compute velocity dependent inertia
%forces for floating frame of reference formulation. Version based on
%minimum amount of matrix transformations. This, less error prone.
%   fbody - structure with information about flexible body

B2 = fbody.B2;
elem_idx = fbody.elem_idx;
q0 = fbody.q0;
n_elem = fbody.n_elem;
elem_len = fbody.elem_len;

% L = [I -A*us A*S]
qn = q0 + B2 * qf;
h = [om; qf_p];

A = Rot(p);
oms = sksym(om);

Qv = zeros(length(q0), 1);
% Loop over elements!
for ii = 1 : n_elem
    qv = @(xi, eta, zeta) Qv_int_L(xi, eta, zeta, elem_idx(:, ii));
    Qv = Qv + volumetric_integral(qv, fbody.mat, fbody.sec, elem_len, 4, 2, 2);
end

function Qv_int = Qv_int_L(xi, eta, zeta, elem_idx)
    S_elem = beam3d_shape_fun(xi, eta, zeta, elem_len);
    
    u = S_elem * qn(elem_idx);
    us = sksym(u);
    
    phi = S_elem * B2(elem_idx, :);
    u_p = phi * qf_p;
    us_p = sksym(u_p);
    
    L = [eye(3), -A * us, A * phi];
    L_p = [-A * (oms * us + us_p), A * oms * phi];
    Qv_int = -L' * (L_p * h);
end

end

