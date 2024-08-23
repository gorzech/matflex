function Fg = beam3d_grav(rho, A, l, g, lumped)
% Compute force due to gravity for 3D beam element based on 
% Shabana and Przemieniecki

if nargin < 5
    lumped = false;
end

m = rho * A * l;

% Integral of the shape matrix
S_int_T = [0.5 * eye(3); zeros(3); 0.5 * eye(3); zeros(3)];
if ~lumped
    l_12 = l / 12;
    S_int_T(5, 3) = -l_12;
    S_int_T(11, 3) = l_12;
    S_int_T(6, 2) = l_12;
    S_int_T(12, 2) = -l_12;
end
    
Fg = S_int_T * g .* m;