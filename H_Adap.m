function [vT, ST, v_paths, s_paths, t_paths, nsteps] = H_Adap(S0, r0, v0, kappa, vbar, gamma, rho_xv, tol, hmax, T, M, keep_paths)
%   此处显示详细说明
v = v0 * ones(M,1);
x = log(S0) * ones(M,1);
t = zeros(M,1); % 每条路当前时刻
nsteps = zeros(M,1); % 每条路的步数

t_paths = {};
v_paths = {};
s_paths = {};
if keep_paths
    v_paths = cell(M,1); 
    t_paths = cell(M,1); 
    s_paths = cell(M,1); 
    for k = 1:M
        t_paths{k} = 0; 
        v_paths{k} = v0; 
        s_paths{k} = S0; 
    end
end

active = true(M,1);
while any(active)
    idx = find(active); % 只处理没到T的路径
    vp = max(v(idx),0);

    h = Adaptive_Step(vp, kappa, vbar, tol, hmax);
    h = min(h, T - t(idx));

    Zv = randn(numel(idx),1); Zx = randn(numel(idx),1);
    dWv = sqrt(h).*Zv;
    dWx = sqrt(h).*(rho_xv*Zv + sqrt(1-rho_xv^2)*Zx);

    v_old = v(idx);
    v(idx) = v_old + kappa*(vbar - vp).*h + gamma*sqrt(vp).*dWv;
    x(idx) = x(idx) + (r0 - 0.5*vp).*h + sqrt(vp).*dWx;  

    t(idx) = t(idx) + h;
    nsteps(idx) = nsteps(idx) + 1;

    if keep_paths
        for j = 1:numel(idx)
            k = idx(j);
            t_paths{k}(end+1) = t(k);
            v_paths{k}(end+1) = v(k);
            s_paths{k}(end+1) = exp(x(k));
        end
    end
    active(idx) = t(idx) < T - 1e-14; % 浮点容差防止死循环
end
ST = exp(x);  vT = v;
end

