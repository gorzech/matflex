function K = beam3d_stiffness(E, G, A, l, Jx, Iy, Iz)
% Compute stiffness matrix K for 3D beam element based on 
% Przemieniecki, 1968, p. 79, simplifications as in Shabana, 2013, p. 295.

k1 = E * A / l;
k2 = 12 * E * Iz / l^3;
k3 = 12 * E * Iy / l^3;
k4 = G * Jx / l;

k_diag = [k1, k2, k3, ...
    k4, 4 * E * Iy / l, 4 * E * Iz / l];

K = diag([k_diag, k_diag]);

K(7, 1) = -k1;
K(1, 7) = -k1;
K(8, 2) = -k2;
K(2, 8) = -k2;
K(9, 3) = -k3; 
K(3, 9) = -k3; 
K(10, 4) = -k4;
K(4, 10) = -k4;

K(11, 5) = 2 * E * Iy / l;
K(5, 11) = K(11, 5);
K(12, 6) = 2 * E * Iz / l;
K(6, 12) = K(12, 6);

k_z_ry = -6 * E * Iy / l^2;
k_y_rz = 6 * E * Iz / l^2;

K(5, 3) = k_z_ry;
K(3, 5) = k_z_ry;
K(6, 2) = k_y_rz;
K(2, 6) = k_y_rz;

K(8, 6) = -k_y_rz;
K(6, 8) = -k_y_rz;
K(9, 5) = -k_z_ry;
K(5, 9) = -k_z_ry;
K(11, 3) = k_z_ry;
K(3, 11) = k_z_ry;
K(12, 2) = k_y_rz;
K(2, 12) = k_y_rz;

K(11, 9) = -k_z_ry;
K(9, 11) = -k_z_ry;
K(12, 8) = -k_y_rz;
K(8, 12) = -k_y_rz;