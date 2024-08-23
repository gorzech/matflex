% Test file for beam3d shape function in order to check it basic properties

% preconditioner - code works and matrix has proper shape
S_pre = beam3d_shape_fun(0, 0, 0, 1);
assert( all(size(S_pre) == [3, 12]) )

%% Check start and end points positions

P1 = [1; 2; 3.5];
P2 = [7; 14; -6];

q = [P1; P1; P2; P2];

L = 1.7;

S1 = beam3d_shape_fun(0, 0, 0, L);
S2 = beam3d_shape_fun(1, 0, 0, L);

% All shape functions should meet value at nodes
assert( norm(P1 - S1 * q) < 1e-14 )
assert( norm(P2 - S2 * q) < 1e-14 )

% Shape functions should always be equal one
assert( norm( sum(S1, 2) - [1; 1; 1] ) < 1e-15 )
assert( norm( sum(S2, 2) - [1; 1; 1] ) < 1e-15 )


%% Check points of undeformed beam

% in beam coordinates they are along x axis. 
P1 = [1; 3; 4];
P2 = [3; 3; 4];
O = zeros(3, 1);

q = [P1; O; P2; O];

L = 2; 

% middle point 
S_mid = beam3d_shape_fun(0.5, 0, 0, L);
assert( norm( S_mid * q - [2; 3; 4]) < 1e-14 )
% 1/3th of the beam at (0.5, 0.5) from centerline
S_13 = beam3d_shape_fun(1/3, 0.5 / L, 0.4 / L, L);
assert( norm( S_13 * q - [1 + 2/3; 3.0; 4.0]) < 1e-14 )

% Shape functions should always be equal one (at least at center-line?)
assert( norm( sum(S_mid, 2) - [1; 1; 1] ) < 1e-15 )

%% Torsion tests - beam should remain straight

P1 = zeros(3, 1);
P2 = [1; 0; 0];
phi1 = 0.2;
phi2 = -0.4;

q = [P1; phi1; 0; 0; P2; phi2; 0; 0];
L = 1;

S_mid = beam3d_shape_fun(0.5, 0, 0, L);
assert( norm( S_mid * q - [0.5; 0; 0]) < 1e-14 )
S_13 = beam3d_shape_fun(1/3, 0.5 / L, 0.4 / L, L);
assert( norm( S_13 * q - [1/3; 0.0; 0.0]) < 1e-14 )

%% Bending in XZ - change angle Y

P1 = zeros(3, 1);
P2 = [1; 0; 0];
L = 1;
% second bending - middle point at centerline
theta1 = 0.05;
theta2 = 0.05;

q = [P1; 0; theta1; 0; P2; 0; theta2; 0];

S_mid = beam3d_shape_fun(0.5, 0, 0, L);
assert( norm( S_mid * q - [0.5; 0; 0]) < 1e-14 )
S_14 = beam3d_shape_fun(1/4, 0.5 / L, 0.4 / L, L);
u = S_14 * q;
assert( abs( u(1) - 0.25 ) < 0.01 && abs( u(2) ) < 1e-14 && u(3) < 0 )
S_34 = beam3d_shape_fun(3/4, 0.5 / L, 0.4 / L, L);
u2 = S_34 * q;
assert( norm( u2 - [0.5 + u(1); 0; -u(3)]) < 1e-14 )

% add some torsion
q(4) = 0.05;
u1 = beam3d_shape_fun(1/4, -0.01 / L, 0.009 / L, L) * q;
u2 = beam3d_shape_fun(1/4, 0.0, 0.009 / L, L) * q;
u3 = beam3d_shape_fun(1/4, 0.01 / L, 0.009 / L, L) * q;
u4 = beam3d_shape_fun(1/4, 0.0 / L, -0.009 / L, L) * q;
% For torsion y coordinate is negative for positive z (indep on y)
assert(u1(2) < 0 && abs(u2(2) - u1(2)) < 1e-14 && abs(u3(2) - u1(2)) < 1e-14)
assert(abs(u4(2) + u1(2)) < 1e-14) % and symmetric
% For torsion z and different y values
assert( abs(u2(3) - u(3)) < 1e-14 )
assert( u2(3) > u1(3) && u2(3) < u3(3) )

%% Bending in XY - change angle Z
P1 = zeros(3, 1);
P2 = [1; 0; 0];
L = 1;
% second bending - middle point at centerline
psi1 = 0.05;
psi2 = 0.05;

q = [P1; 0; 0; psi1; P2; 0; 0; psi2];

S_mid = beam3d_shape_fun(0.5, 0, 0, L);
assert( norm( S_mid * q - [0.5; 0; 0]) < 1e-14 )
S_14 = beam3d_shape_fun(1/4, 0.5 / L, 0.4 / L, L);
u = S_14 * q;
assert( abs( u(1) - 0.25 ) < 0.01 && abs( u(3) ) < 1e-14 && u(2) > 0 )
S_34 = beam3d_shape_fun(3/4, 0.5 / L, 0.4 / L, L);
u2 = S_34 * q;
assert( norm( u2 - [0.5 + u(1); -u(2); 0]) < 1e-14 )

% add some torsion
q(4) = -0.05;
u1 = beam3d_shape_fun(1/4, -0.01 / L, 0.009 / L, L) * q;
u2 = beam3d_shape_fun(1/4, -0.01 / L, 0.0, L) * q;
% For torsion z coordinate is positive for negative y (indep on z)
assert(u1(3) > 0 && abs(u2(3) - u1(3)) < 1e-14)
% For torsion y and different z values
assert( abs(u2(2) - u(2)) < 1e-14 )
assert( u2(2) < u1(2) )