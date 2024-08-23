%% Basic test to evaluate approxJacobian

x0 = [1.23, 6.15, 9.0]';

jac = approxJacobian(@vector_test_function, x0);

jac_exact = exact_jacobian(x0);

jac_norm = norm(jac_exact - jac);
assert(jac_norm < 2.4e-6);


function q = vector_test_function(x) 
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    
    q = [0.0
        146.0
    	x1 + 4.50 * x2 * x2 + sqrt(x3) / 12.0
        13.0 * sin(4.0 * x1) + cos(3.0 * x2)
        1.0 / (x1 + x2 + x3)
        exp(2.0 * x1) + log(x2)];
end

function m = exact_jacobian(x) 
    row_4 = -1.0 / (x(1) + x(2) + x(3)) ^ 2.0;
    m = [0.0, 0.0, 0.0
        0.0, 0.0, 0.0 
        1.0, 9.0 * x(2), 1.0 / 24.0 / sqrt(x(3))
        52.0 * cos(4.0 * x(1)), -3.0 * sin(3.0 * x(2)), 0.0
        row_4, row_4, row_4
        2.0 * exp(2.0 * x(1)), 1.0 / x(2), 0.0];
end