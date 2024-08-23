function S = beam3d_shape_fun(xi, eta, zeta, L)
%BEAM3D_SHAPE_FUN Compute shape function of 3d beam element
%   Compute shape function for a 3d beam element with 12 nodal coordinates
%   based on book (Shabana, 2013) Dynamics of multibody systems
%   xi, eta, zeta are dimentionless spacial parameters
%   l is element length

xi2 = xi .^ 2;
xi3 = xi2 .* xi;
w1 = 1 - xi;
w2 = 6 * xi * w1;
w3 = (1 - 4 * xi + 3 * xi2) * L;
w4 = (3 .* xi2 - 2 .* xi) * L;
w5 = 1 - 3 * xi2 + 2 * xi3;
w6 = (xi - 2 * xi2 + xi3) * L;

S = [w1,    w2 * eta,   w2 * zeta,  0,              w3 * zeta,  -w3 * eta, ...
        xi, -w2 * eta,  -w2 * zeta, 0,              w4 * zeta,      -w4 * eta
    0,      w5,         0,          -w1 * L * zeta, 0,          w6, ...
        0,  1 - w5,     0,          -L * xi * zeta, 0,              -L * xi2 * w1
    0,      0,          w5,         w1 * L * eta,  -w6,        0, ...
        0,  0,          1 - w5,     L * xi * eta,  L * xi2 * w1,   0];

end

