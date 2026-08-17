function [FixedTime, AdaptiveTime, AdaptiveSteps] = ...
    measure_cost(P, T_now, N_now, hmax_now, tol_now, R)

FixedRuns = nan(R,1);
AdaptiveRuns = nan(R,1);
AdaptiveStepRuns = nan(R,1);

for j = 1:R

    % Fixed-step Euler
    rng(P.seed + 1000*j);

    tic;
    [ST_F, I_F, ~] = HCIR_Fixed( ...
        P.S0, P.r0, P.lambda, P.theta, P.eta, ...
        P.v0, P.kappa, P.vbar, P.gamma, ...
        P.rho_xv, P.rho_xr, N_now, P.M, T_now, false);

    payoffF = exp(-I_F(:)) .* max(ST_F(:)-P.K, 0);
    mean(payoffF); %#ok<VUNUS>
    FixedRuns(j) = toc;

    % Adaptive-step Euler
    rng(P.seed + 2000*j);

    tic;
    [ST_V, I_V, ~, ~, ~, ns] = HCIR_Adaptive( ...
        P.r0, P.lambda, P.theta, P.eta, ...
        P.v0, P.kappa, P.vbar, P.gamma, ...
        tol_now, hmax_now, P.S0, ...
        P.rho_xv, P.rho_xr, P.M, T_now, false);

    payoffV = exp(-I_V(:)) .* max(ST_V(:)-P.K, 0);
    mean(payoffV); %#ok<VUNUS>
    AdaptiveRuns(j) = toc;
    AdaptiveStepRuns(j) = mean(ns);
end

% Median is less affected by occasional background CPU activity
FixedTime = median(FixedRuns);
AdaptiveTime = median(AdaptiveRuns);
AdaptiveSteps = mean(AdaptiveStepRuns);
end

