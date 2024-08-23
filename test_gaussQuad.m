% Verify if basic integration for gaussQuad works as expected

%% Integrate constant 

value = 1.17;
low = -7.12;
hi = 3.14;

test = @(~) value;

result = gaussQuad(test, low, hi, 3);

expected = value * (hi - low);

assert( abs(result - expected) < 4e-15 )

%% Test for linear function 

% Taken from the book 
% Svein Linge and Hans Petter Langtangen, Programming for Computations – MATLAB/Octave

f = @(x) 6*x - 4;
F = @(x) 3*x^2 - 4*x; % Anti-derivative
a = 1.2; b = 4.4;
expected = F(b) - F(a);
tol = 2e-14;
for n = [2 20 21]
    computed = gaussQuad(f, a, b, n);
    error = abs(expected - computed);
    assert(error < tol, 'n=%d, err=%g', n, error);
end