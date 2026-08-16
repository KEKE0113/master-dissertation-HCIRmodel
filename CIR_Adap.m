function [rT, r_paths, t_paths, nsteps] = CIR_Adap(r0, lambda, theta, eta, tol, hmax, T, M, keep_paths)
% CIR_ADAP Simulates the CIR short-rate process using an adaptive full truncation Euler scheme.
% Model: dr(t) = lambda * (theta - r(t)) dt + eta * sqrt(r(t)) dW(t)

% Initialise the short rate and time
r = r0 * ones(M, 1);
t = zeros(M, 1);
nsteps = zeros(M, 1);

% Initialise path storage if required
if keep_paths
    r_paths = cell(M, 1);
    t_paths = cell(M, 1);

    for m = 1:M
        r_paths{m} = r0;
        t_paths{m} = 0;
    end
else
    r_paths = {};
    t_paths = {};
end

% Identify paths that have not yet reached maturity
active = true(M,1);

while any(active)

    idx = find(active); 

    % Determine the adaptive time step for each active path
    h = Adaptive_Step(r(idx), lambda, theta, tol, hmax);

    % Ensure that no step moves beyond maturity
    h = min(h, T - t(idx));

    dWr = sqrt(h) .* randn(size(h));

    r_old = r(idx);
    r_pos = max(r_old,0);

    r(idx) = r_old + lambda * (theta - r_pos) .* h + eta * sqrt(r_pos) .* dWr;

    t(idx) = t(idx) + h;
    nsteps(idx) = nsteps(idx) + 1;

    if keep_paths
        for j = 1:numel(idx)
            m = idx(j);

            t_paths{m}(end + 1) = t(m);
            r_paths{m}(end + 1) = r(m);
        end
    end

    % Remove paths that have reached maturity
    active(idx) = t(idx) < T - 1e-14;
end

rT = r;

end

