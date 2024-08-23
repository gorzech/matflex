% test_gauleg - verify if gauleg produces reasonable results
% As I believe I got this file from the Internet

% Compare with exact results provided by wikipedia:
% https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature

%% Test for 1 point

% For default interval
[x, w] = gauleg(-1, 1, 1);
assert( isscalar(x) && isscalar(w) )
assert( norm(x - 0) < 1e-15 )
assert( norm(w - 2) < 1e-15 )

% for custom interval
[x, w] = gauleg(5, 6, 1);
assert( isscalar(x) && isscalar(w) )
assert( norm(x - 5.5) < 1e-15 )
assert( norm(w - 1) < 1e-15 )

%% Test for four points

w_1_4 = (18 - sqrt(30)) / 36;
w_2_3 = (18 + sqrt(30)) / 36;

x_1_4_abs = sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5));
x_2_3_abs = sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5));
[x, w] = gauleg(-1, 1, 4);

x_expected = [-x_1_4_abs, -x_2_3_abs, x_2_3_abs, x_1_4_abs];
w_expected = [w_1_4, w_2_3, w_2_3, w_1_4];

assert( norm(x(:) - x_expected(:)) < 1e-15 )
assert( norm(w(:) - w_expected(:)) < 1e-15 )