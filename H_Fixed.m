function [ST, vT, spaths, vpaths] = H_Fixed(S0, r0, v0, kappa, vbar, gamma, T, N, M, keep_paths)
%H_FIXED 用HESTON Model模拟在固定步长下的v(t)
%   输入: kappa:均值回归速度；vbar:均值长期水平；gamma:v(t)的波动率
%   输入：T：时间； N：时间步数； M：模拟次数

dt = T/N; % 时间步长
x = log(S0) * ones(M,1); % log price
v = v0 * ones(M,1); % variance(总共有M个模拟结果)

if keep_paths % 是否保存路径
    spaths = zeros(M, N+1); spaths(:,1) = S0;
    vpaths = zeros(M, N+1); vpaths(:,1) = v;
end

for i = 1:N % 时间维度（第i次迭代）
    v_pos = max(0,v);

    Zv = randn(M,1); Zx = randn(M,1);
    dWv = sqrt(dt) * Zv;
    dWx = sqrt(dt) * Zx;
    
    v = v + kappa*(vbar - v_pos)*dt + gamma*sqrt(v_pos).*dWv;
    x = x + (r0 - 0.5*v_pos)*dt + sqrt(v_pos).*dWx;
    
    if keep_paths
       vpaths(:,i+1) = v;
       spaths(:,i+1) = exp(x);
    end % 保存每一次路径
end

ST = exp(x);
vT = v;
end


