function [V, F] = normalized_eig(K, M, n_modes)

% In general try to return sorted frequencies and M-normalized Vectors
if issparse(K)
    K = full(K);
end
if issparse(M)
    M = full(M);
end
if nargin < 3
    n_modes = size(K, 1);
end

[eigV, D] = eig(K, M, 'vector'); 

% If e.g. M is singular, we get NaN and Inf values. Then 
numbers_D = find(isfinite(D)); % return indices

n_modes = min(n_modes, length(numbers_D));

[d, idx] = sort(D(numbers_D));

if nargin > 1
    F = sqrt(abs(d(1 : n_modes))) / 2 / pi;
end

% get target V
V = eigV(:, numbers_D(idx(1 : n_modes)));

% Now need to M - normalize them
V = V * diag( sqrt( 1 ./ diag(V'*(M*V)) ) );