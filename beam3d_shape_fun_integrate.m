function S = beam3d_shape_fun_integrate(rho, A, L)
%BEAM3D_SHAPE_FUN Compute integrated shape function of 3d beam element

w1 = L / 12;
S = [0.5,    0,   0,  0,              0,  0, ...
        0.5, 0,  0, 0,              0,      0
    0,      0.5,         0,          0, 0,          w1, ...
        0,  0.5,     0,          0, 0,              -w1
    0,      0,          0.5,         0,  -w1,        0, ...
        0,  0,          0.5,     0,  w1,   0] .* (rho * A * L);

end

