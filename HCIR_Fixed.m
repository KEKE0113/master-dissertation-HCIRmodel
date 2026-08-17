function [ST, I, paths] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N, M, T, keep_paths)
%HESTON_CIR_FIXED 模拟在hybrid model下固定步长的S(t)
%   此处采用对数形式模拟 x = ln S(t)

dt = T/N; % 时间步长
x = log(S0) * ones(M,1); % 总共有M个模拟结果
r = r0 * ones(M,1);
v = v0 * ones(M,1);
I  = zeros(M,1); 

c3 = sqrt(max(1 - rho_xv^2 - rho_xr^2, 0)); % Cholesky, (x,v,r) 排序

if keep_paths % 是否保存路径
    paths = zeros(M, N+1);
    paths(:,1) = exp(x);
else
    paths = [];
end

for i = 1:N % 时间维度（第i次迭代）
    r_old = r; v_old = v; x_old = x;
    rp = max(r_old, 0); vp = max(v_old, 0);

    Z = randn(M,3); Zx = Z(:,1); Zv = Z(:,2); Zr = Z(:,3);
    a = sqrt(1 - rho_xv^2);
    dWx = sqrt(dt) * Zx;
    dWv = sqrt(dt) * (rho_xv * Zx + a * Zv);
    dWr = sqrt(dt) * (rho_xr * Zx - (rho_xv * rho_xr / a) * Zv ...
    + sqrt((1 - rho_xv^2 - rho_xr^2) / (1 - rho_xv^2)) * Zr);

    x = x_old + (rp-0.5*vp)*dt + sqrt(vp).*dWx;
    r = r_old + lambda*(theta - rp)*dt + eta*sqrt(rp).*dWr;
    v = v_old + kappa*(vbar - vp)*dt + gamma*sqrt(vp).*dWv;

    I = I + 0.5*(rp + max(r,0))*dt;
    if keep_paths, paths(:,i+1) = exp(x); end % 保存每一次路径
end

ST = exp(x);
end

