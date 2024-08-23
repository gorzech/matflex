function int_rho_dV = volumetric_integral(intfun, mat, sec, len, rankx, ranky, rankz)
%VOLUMETRIC_INTEGRAL Helper function that uses triple gauss quad
%   Detailed explanation goes here
half_eta = sec.thky / 2 / len;
half_zeta = sec.thkz / 2 / len;
int_rho_dV = tripleGaussQuad(intfun, 0, 1, -half_eta, half_eta, -half_zeta, half_zeta, ...
    rankx, ranky, rankz);
int_rho_dV = int_rho_dV .* (mat.rho * len ^ 3);
end

