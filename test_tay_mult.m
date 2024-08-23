% Test script to verify if tay_mult works as expected (up to order1)

% Generate test data

nrow = 3; 
ncol = 4;
nq = 7;

M0 = rand(nrow, ncol);

M1_1 = rand(nrow, nq);
M1_2 = rand(nrow, nq);
M1_3 = rand(nrow, nq);
M1_4 = rand(nrow, nq);

q = rand(7, 1);

assert(isequal(M0, u_tay_mult(M0)))
assert(isequal(M0, u_tay_mult(M0, q)))

%% Test more exotic types

m.m0 = M0;
m.m1 = [];

assert(isequal(M0, u_tay_mult(m)))

m.m1 = cat(3, M1_1, M1_2, M1_3, M1_4);

expected = M0 + [M1_1 * q, M1_2 * q, M1_3 * q, M1_4 * q];
assert(norm(expected - u_tay_mult(m, q)) < 1e-15)

m.m0 = [];
expected = expected - M0;
assert(norm(expected - u_tay_mult(m, q)) < 1e-15)