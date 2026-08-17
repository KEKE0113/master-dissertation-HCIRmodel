%% Parameters
S0 = 100;
K =  95;
r0 = 0.03;
lambda = 0.5;
theta = 0.04; 
eta = 0.15;
v0 = 0.04;
kappa = 0.3;
vbar = 0.04;
gamma = 0.15;
rho_xv = -0.7;
rho_xr = 0.6;
N = 200; % Number of fixed steps
M = 100000; % Number of simulations
T = 1; 
tol = 0.00001;
hmax = T/200; % 0.005
seed = 20260804;

%% Fixed Steps
rng(seed);
[ST_F, I_F, S_paths_F] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N, M, T, true);

% stock price figure
t = linspace(0,T,N+1);
num_paths = 20;
idx_F = randperm(M, num_paths); % 随机选10条路径

figure;
plot(t,S_paths_F(idx_F,:)','LineWidth',1)
xlabel('Time')
ylabel('Stock Price'); 
title('Stock Price Paths under Heston-CIR model with fixed step size')
set(gca, 'FontSize', 12 )
grid on


%% Adaptive 
rng(seed);
[ST_V, I_V, t_paths_V, S_paths_V, h_paths_V, nsteps_V] = HCIR_Adaptive(r0, lambda, theta, eta, v0, kappa, vbar, gamma, tol, hmax, S0, rho_xv, rho_xr, M, T, true);

% stock price figure
num_paths = 20;
idx_V = randperm(M, num_paths); % 随机选10条路径

figure; hold on;
for k = idx_V
    plot(t_paths_V{k}, S_paths_V{k}, 'LineWidth', 1);
end
xlim([0 T]);
xlabel('Time'); 
ylabel('Stock Price'); ylim([60 180]);
title('Stock Price Paths under Heston-CIR model with variable step size');
grid on; hold off;

% Average step size and the number of steps
avg_steps = mean(nsteps_V); % 平均的步数


%% option price figure (option price with K change)

K_test = [65, 80, 95, 110, 125];
C_F = zeros(1,length(K_test));
C_V = zeros(1,length(K_test));
C_b = zeros(1,length(K_test));
SE_F = zeros(1,length(K_test));
SE_V = zeros(1,length(K_test));

for i = 1:length(K_test)
    payoff_V = exp(-I_V) .* max(ST_V-K_test(i),0);
    payoff_F = exp(-I_F) .* max(ST_F-K_test(i),0);
    C_V(i) = mean(payoff_V);
    C_F(i) = mean(payoff_F);
    C_b(i) = HCIR_ChF_Call(S0, K_test(i), T, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr);

    SE_F(i) = std(payoff_F) / sqrt(M);
    SE_V(i) = std(payoff_V) / sqrt(M);
    dev_f(i) = C_F(i) - C_b(i); 
    dev_a(i) = C_V(i) - C_b(i);
    CI_fl(i) = C_F(i)-1.96*SE_F(i); CI_fh(i) = C_F(i)+1.96*SE_F(i); 
    CI_al(i) = C_V(i)-1.96*SE_V(i); CI_ah(i) = C_V(i)+1.96*SE_V(i);
end

figure; hold on;
plot(K_test,C_F,'o-','LineWidth',1)
plot(K_test,C_V,'o-','LineWidth',1)
plot(K_test,C_b,'o-','LineWidth',1)
legend('Fixed step', 'Adaptive', 'Benchmark');
xlabel('Strike Price')
ylabel('Call Option Price')
title('Call option price with different strike price')
set(gca, 'FontSize', 12 )
grid on; hold off;

% Error vs K
figure; hold on;
plot(K_test,abs(dev_f),'LineWidth',1)
plot(K_test,abs(dev_a),'LineWidth',1)
legend('Fixed step', 'Adaptive step');
xlabel('Strike')
ylabel('Deviation')
title('Deviation from the reference price with different strike price')
set(gca, 'FontSize', 12 )
grid on; hold off;

fprintf('\n%-4s %9s | %9s %7s %9s %6s %-4s | %9s %7s %9s %6s %-4s\n',...
   'K','C_ref','C_fix','SE','dev','CI','C_adp','SE','dev','CI');
for i = 1:length(K_test)
    fprintf('\n%-4d %9.4f | %9.4f %7.4f %+9.4f [%8.4f,%8.4f] | %9.4f %7.4f %+9.4f [%8.4f,%8.4f]\n',...
        K_test(i), C_b(i),...
        C_F(i),SE_F(i),dev_f(i),CI_fl(i),CI_fh(i),...
        C_V(i),SE_V(i),dev_a(i),CI_al(i),CI_ah(i));
end



%% Deviation vs M (K=95) 模拟次数
C_b1 = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr);
M_test = [1000, 5000, 10000, 50000, 100000]; % 模拟的次数
N_test = 30; % 每个M重复N次

err_meanF = zeros(size(M_test)); err_meanV= zeros(size(M_test)); % 平均绝对误差
SE_meanF = zeros(size(M_test)); SE_meanV = zeros(size(M_test)); % 平均标准误
cover_F = zeros(size(M_test)); cover_V = zeros(size(M_test));% 95% CI覆盖benchmark的比例

for j = 1:length(M_test) 
    n = M_test(j); 
    err_V = zeros(1,N_test); err_F = zeros(1,N_test);
    SEs_V = zeros(1,N_test); SEs_F = zeros(1,N_test);
    hit_V = zeros(1,N_test); hit_F = zeros(1,N_test);
    for m = 1:N_test 
        [ST_V, I_V, t_paths, S_paths, h_paths, nsteps] = HCIR_Adaptive(r0, lambda, theta, eta, v0, kappa, vbar, gamma, tol, hmax, S0, rho_xv, rho_xr, n, T, false);
        [ST_F, I_F, paths] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N, n, T, false);
        
        payoff_V = exp(-I_V) .* max(ST_V-K,0); 
        payoff_F = exp(-I_F) .* max(ST_F-K,0);

        C_V = mean(payoff_V); 
        C_F = mean(payoff_F);

        SE_V = std(payoff_V)/sqrt(n); 
        SE_F = std(payoff_F)/sqrt(n);

        err_V(m) = C_V - C_b1; 
        err_F(m) = C_F - C_b1;

        SEs_V(m) = SE_V; 
        SEs_F(m) = SE_F; 

        hit_V(m) = (C_b1 >= C_V-1.96*SE_V) && (C_b1 <= C_V+1.96*SE_V);
        hit_F(m) = (C_b1 >= C_F-1.96*SE_F) && (C_b1 <= C_F+1.96*SE_F);
    end
    err_meanF(j) = mean(err_F); 
    err_meanV(j) = mean(err_V);

    SE_meanF(j) = mean(SEs_F); 
    SE_meanV(j) = mean(SEs_V);

    cover_F(j) = mean(hit_F); 
    cover_V(j) = mean(hit_V);
    fprintf('N=%7d | err_F=%.4f | SE_F=%.4f | err/SE_F = %.4f\n', ...
            n, err_meanF(j), SE_meanF(j), abs(err_meanF(j)/SE_meanF(j)));
    fprintf('N=%7d | err_V=%.4f | SE_V=%.4f | err/SE_V = %.4f\n', ...
            n, err_meanV(j), SE_meanV(j), abs(err_meanV(j)/SE_meanV(j)));
end

% Plot figures
figure; hold on;
loglog(M_test, SE_meanF, 's-','LineWidth',1.2);
loglog(M_test, SE_meanV, 'd-','LineWidth',1.2);
ref = SE_meanF(1)*sqrt(M_test(1))./sqrt(M_test);   % 1/√M 参考线
loglog(M_test, ref, 'k:','LineWidth',1);
set(gca,'XScale','log','YScale','log','FontSize',12);
xlabel('Number of Simulations'); ylabel('Standard error (sampling noise)');
legend('Fixed','Adaptive','Slope -1/2','Location','southwest');
title('Sampling noise vs M'); grid on; hold off;


%% Computational Time

% 平均总步数 vs bias
N_fixed = [100, 400, 800, 1600, 3200];
Tol1 = [0.005, 0.001, 0.0005, 0.0002, 0.0001];
hmax = T/10;
R = 20;

C_b = zeros(1,length(N_fixed));
C_F = zeros(1,length(N_fixed));
Bias_F = zeros(1,length(N_fixed));
for i = 1:length(N_fixed) % 求出对应N_fixed步数下C_F的bias
    for r = 1:R
        rng(1000*r)
        C_b(i) = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr);
        [ST, I, paths] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N_fixed(i), M, T, false);
        payoff_F = exp(-I) .* max(ST-K,0);
        C_F(i) = mean(payoff_F);
    end
    Bias_F(i) = C_F(i)-C_b(i);
end

C_V = zeros(1,length(N_fixed));
Steps = zeros(1,length(N_fixed));
Bias_V = zeros(1,length(N_fixed));
for i = 1:length(Tol1) % 求出Tol对应步数下C_V的bias
    for r = 1:R
        rng(1000*r);
        [ST, I, ~, ~, ~, nsteps] = HCIR_Adaptive(r0, lambda, theta, eta, v0, kappa, vbar,gamma, Tol1(i), hmax, ...
        S0, rho_xv, rho_xr, M, T, false);
        Steps(i) = mean(nsteps);
        payoff_V = exp(-I) .* max(ST-K,0);
        C_V(i) = mean(payoff_V);
        Bias_V(i) = C_V(i)-C_b(i);
    end
end

figure; hold on;
plot(N_fixed, Bias_F, 'LineWidth',1);
plot(Steps, Bias_V, 'LineWidth',1);
xlabel('Number of Steps'); ylabel('Bias');
legend('Fixed step', 'Variable step');
title('Bias vs Steps')
set(gca, 'FontSize', 12 )
grid on; hold off;

%% 扫步长（没有频繁触0）
N_list = [50 200 800 1600 2000];  % 几何等比
h_list = T ./ N_list;
bias_F = zeros(size(N_list));
biasSE_F = zeros(size(N_list));
R = 20;

C_ref1 = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr);
for j = 1:numel(N_list)
    bia_F = zeros(1,R);
    for r = 1:R
        rng(seed+r);                              % CRN across N
        [ST,I,~] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N_list(j), M, T, false);
        bia_F(r) = mean(exp(-I).* max(ST-K,0)) - C_ref1;
    end
    bias_F(j) = mean(bia_F);
    biasSE_F(j) = std(bia_F)/sqrt(R); 
    fprintf('N=%5d, h=%.4f | bias=%+.4e ± %.1e\n', ...
            N_list(j), T/N_list(j), bias_F(j), biasSE_F(j));
end
figure;
loglog(h_list, abs(bias_F), 'o-'); hold on
loglog(h_list, abs(bias_F(1))*(h_list/h_list(1)), 'k--')  % slope-1 参考
xlabel('Step size h'); ylabel('|C(h) - C_{ref}|')
legend('Measured bias','Slope 1 (weak Euler)','Location','northwest')
title('Discretization convergence under fixed step'); grid on
set(gca,'FontSize',12)
% 拟合斜率
p = polyfit(log(h_list), log(abs(bias_F)), 1);
fprintf('fixed weak order ≈ %.2f\n', p(1));   % 期望 ≈1

%% 扫步长（频繁触0）
gamma=0.3;
N_list = [50 100 200 400 800];  % 几何等比
h_list = T ./ N_list;
bias_F = zeros(size(N_list));
biasSE_F = zeros(size(N_list));
R = 20;

C_ref1 = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr);
for j = 1:numel(N_list)
    bia_F = zeros(1,R);
    for r = 1:R
        rng(seed+r);                              % CRN across N
        [ST,I,~] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, N_list(j), M, T, false);
        bia_F(r) = mean(exp(-I).* max(ST-K,0)) - C_ref1;
    end
    bias_F(j) = mean(bia_F);
    biasSE_F(j) = std(bia_F)/sqrt(R); 
    fprintf('N=%5d, h=%.4f | bias=%+.4e ± %.1e\n', ...
            N_list(j), T/N_list(j), bias_F(j), biasSE_F(j));
end
figure;
loglog(h_list, abs(bias_F), 'o-'); hold on
loglog(h_list, abs(bias_F(1))*(h_list/h_list(1)), 'k--')  % slope-1 参考
xlabel('Step size h'); ylabel('|C(h) - C_{ref}|')
legend('Measured bias','Slope 1 (weak Euler)','Location','northwest')
title('Discretization convergence under fixed step(gamma = 0.3)'); grid on
set(gca,'FontSize',12)
% 拟合斜率
p = polyfit(log(h_list), log(abs(bias_F)), 1);
fprintf('fixed weak order ≈ %.2f\n', p(1));   % 期望 ≈1



%%
T1 = [1, 2, 5, 10, 20]; % 时间长度与效率的关系
C_b = zeros(1,length(T1));
Err_F1 = zeros(1,length(T1)); Err_V1 = zeros(1,length(T1));
fprintf('T = %.4f',T1);
for i = 1:length(T1)
    t = T1(i);
    C_b(i) = HCIR_ChF_Call(S0, K, t, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr);
    [t_F, C_F, t_V, C_V, avg_VSteps] = HCIR_Eff(K, r0, lambda, theta, eta, ...
    v0, kappa, vbar, gamma, ...
    tol, hmax, ...
    S0, rho_xv, rho_xr, M, N1, t);
    Err_F1(i) = abs(C_F-C_b(i));
    Err_V1(i) = abs(C_V-C_b(i));
    fprintf('Fsteps=%7d | time_Fixed=%.4f | time_Adap=%.4f | Vsteps=%.4f\n' , ...
            N1, t_F, t_V, avg_VSteps);
end

[C_b2, P0T] = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr);
T = 5;
N1 = 1000;
Tol = [0.01, 0.005, 0.001, 0.0005, 0.0001];
Err_V2 = zeros(1,length(Tol));
fprintf('Tolerance = %.4f',Tol);
for i = 1:length(Tol)
    tol = Tol(i);
    [t_F, C_F, t_V, C_V, avg_VSteps] = HCIR_Eff(K, r0, lambda, theta, eta, ...
    v0, kappa, vbar, gamma, ...
    tol, hmax, ...
    S0, rho_xv, rho_xr, M, N1, T);
    fprintf('Fsteps=%7d | time_Fixed=%.4f | time_Adap=%.4f | Vsteps=%.4f\n' , ...
            N1, t_F, t_V, avg_VSteps);
    Err_V2(i) = abs(C_V - C_b2);
end




