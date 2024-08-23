function [sec, mat, beam_len, n_elem] = th_beam_fem_props()
% 2 m flexible beam - for different test purposes

% Define material mat and section properties sec
mat.rho = 7801;
beam_len = 2;
% Thin-walled cross section
sec.A = 0.006144;
sec.Iy = 3.7814e-05;
sec.Iz = 3.7814e-05;
sec.Jx = 4.9836032e-05; % sec.Iy + sec.Iz;

mat.E = 2e11;
nu = 0.3;
mat.G = mat.E / 2 / (1 + nu);

n_elem = 40;

