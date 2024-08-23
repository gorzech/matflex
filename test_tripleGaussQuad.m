% simple tests to check in tripleGaussQuad works

%% Integrate constant over volume

value = -14.2546;

low1 = -0.2352;
hi1 = -0.234;
low2 = -102;
hi2 = 9;
low3 = 1;
hi3 = 2;

f = @(~, ~, ~) value;

result = tripleGaussQuad(f, low1, hi1, low2, hi2, low3, hi3, 1, 3, 2);

expected = value * (hi1 - low1) * (hi2 - low2) * (hi3 - low3);

assert( abs(result - expected) < 1e-15 )

%% Another pretty simple example

expected = 31 / 12;

f = @(x, y, z) x * x + y * z;

result = tripleGaussQuad(f, 0, 1, 1, 2, 1, 2,  2, 1, 1);

assert( abs(result - expected) < 1e-15 )

%% Once again - non polynomial function

expected = 180;
f = @(x, y, z) x * x * y * (cos(pi * z) + 2);

result = tripleGaussQuad(f, 0, 3, 0, 2, 0, 5,  2, 1, 3);
error = abs(expected - result);
assert(error < 1e-14, 'err=%g', error);