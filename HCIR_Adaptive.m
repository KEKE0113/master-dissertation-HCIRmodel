function [ST, I, t_paths, S_paths, h_paths, nsteps] = HCIR_Adaptive(r0, lambda, theta, eta, ...
    v0, kappa, vbar, gamma, ...
    tol, hmax, ...
    S0, rho_xv, rho_xr, M, T, keep_paths)

%HCIR_ADAPTIVE 此处显示有关此函数的摘要
%   输入CIR: r0, lambda, theta, eta
%   输入Heston: v0, kappa, vbar, gamma

x = log(S0) * ones(M,1); % 总共有M个模拟结果
r = r0 * ones(M,1);
v = v0 * ones(M,1);
I = zeros(M,1); 
t = zeros(M,1); % 每条路当前时刻
nsteps = zeros(M,1); % 每条路的步数
c3 = sqrt(max(1 - rho_xv^2 - rho_xr^2, 0));

t_paths = {}; 
h_paths = {};
S_paths = {}; 
if keep_paths
    t_paths = cell(M,1); 
    S_paths = cell(M,1); 
    for k = 1:M
        t_paths{k} = 0;  
        h_paths{k} = [];
        S_paths{k} = S0; 
    end
end

active = true(M,1);
while any(active)
    idx = find(active); % 只处理没到T的路径

    h_CIR = Adaptive_Step(r(idx), lambda, theta, tol, hmax);
    h_Hes = Adaptive_Step(v(idx), kappa, vbar,tol, hmax);
    h = min(h_CIR, h_Hes);
    h = min(h, T - t(idx));  

    r_old = r(idx); v_old = v(idx);
    rp = max(r_old,0);  vp = max(v_old,0);
    
    Z = randn(numel(idx), 3); Zx = Z(:,1); Zv = Z(:,2); Zr = Z(:,3);
    a = sqrt(1 - rho_xv^2);
    dWx = sqrt(h) .* Zx;
    dWv = sqrt(h) .* (rho_xv * Zx + a * Zv);
    dWr = sqrt(h) .* (rho_xr * Zx - (rho_xv * rho_xr / a) * Zv ...
    + sqrt((1 - rho_xv^2 - rho_xr^2) / (1 - rho_xv^2)) * Zr);

    x(idx) = x(idx) + (rp - 0.5*vp).*h + sqrt(vp).*dWx;
    r(idx) = r_old + lambda*(theta - rp).*h + eta*sqrt(rp).*dWr;   % full truncation
    v(idx) = v_old + kappa*(vbar  - vp).*h + gamma*sqrt(vp).*dWv;
    I(idx) = I(idx) + 0.5*(rp + max(r(idx),0)).*h;

    t(idx) = t(idx) + h;
    nsteps(idx) = nsteps(idx) + 1;

    if keep_paths
        for j = 1:numel(idx)
            k = idx(j);
            t_paths{k}(end+1) = t(k);
            h_paths{k}(end+1) = h(j);
            S_paths{k}(end+1) = exp(x(k));
        end
    end

    active(idx) = t(idx) < T - 1e-14; % 浮点容差防止死循环
end
ST = exp(x);
end
