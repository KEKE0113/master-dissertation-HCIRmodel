function [rT, paths] = CIR_Fixed(r0, lambda, theta, eta, T, N, M, keep_paths)
% CIR_FIXED Simulates the CIR short-rate process using a fixed-step full truncation Euler scheme.
% Model: dr(t) = lambda * (theta - r(t)) dt + eta * sqrt(r(t)) dW(t)

dt = T / N;

% Initialise the short rate
r = r0 * ones(M, 1);

% Preallocate path storage if required
if keep_paths
    paths = zeros(M, N+1);
    paths(:,1) = r;
else
    paths = [];
end

% Fixed-step Euler simulation
for n = 1:N
    r_pos = max(r,0);

    dW = sqrt(dt) * randn(M,1);

    r = r + lambda * (theta - r_pos)*dt + eta * sqrt(r_pos).* dW;
    
    if keep_paths
        paths(:,n+1) = r; 
    end 
end

rT = r;

end
