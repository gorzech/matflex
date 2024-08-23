function C = approxJacobian(fun, er, tol)
%COMPUTEJAKOBIAN compute jakobian based on function and vector of variables
ner = size(er,1);

if nargin < 3
    delt = sqrt(eps) * abs(er); % Approximately sqrt of machine precision
    delt(delt == 0) = sqrt(eps);
else
    if isscalar(tol)
        delt = ones(ner, 1) .* tol;
    else
        delt = tol;
    end
end

f0 = fun(er);
C = zeros(length(f0), ner);
for wi = 1 : ner
    esafe = er(wi);
    er(wi) = esafe + delt(wi);
    f = fun(er);
    C(:, wi) = (f - f0) ./ delt(wi);
    er(wi) = esafe;
end