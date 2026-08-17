function [t_F, C_F, t_V, C_V, avg_VSteps] = HCIR_Eff(K, r0, lambda, theta, eta, ...
    v0, kappa, vbar, gamma, ...
    tol, hmax, ...
    S0, rho_xv, rho_xr, M, N, T)
%HCIR_TIMEF 此处显示有关此函数的摘要
%   此处显示详细说

tic;
[ST_F, I_F] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N, M, T, false);
t_F = toc;
payoff_F = exp(-I_F) .* max(ST_F-K,0);
C_F = mean(payoff_F);


tic;
[ST_V, I_V, ~, ~, ~, nsteps] = HCIR_Adaptive(r0, lambda, theta, eta,v0, kappa, vbar, gamma,tol, hmax, S0, rho_xv, rho_xr, M, T, false);
t_V = toc;
avg_VSteps = mean(nsteps);
payoff_V = exp(-I_V) .* max(ST_V-K,0);
C_V = mean(payoff_V);

end

