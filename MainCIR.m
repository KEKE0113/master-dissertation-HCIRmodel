%% Parameters
S0 = 100;
K =  95;
r0 = 0.04;
lambda = 0.5;
theta = 0.04; 
eta = 0.15;
N = 200; % Number of fixed steps
M = 100000; % Number of simulations
T = 1; 
tol = 0.00001;
hmax = T/200;
seed = 20260731;

%% Path
rng(seed);  [rF, PathsF] = CIR_Fixed(r0, lambda, theta, eta, T, N, M, true);
rng(seed);  [rA, PathsA, tPaths, nsteps] = CIR_Adap(r0, lambda, theta, eta, tol, hmax, T, M, true);

t = linspace(0,T,N+1);
num_paths = 20;
idx = randperm(M, num_paths); % 随机选10条路径

figure; hold on;
plot(t,PathsF(idx,:)','LineWidth',1)
xlabel('Time')
ylabel('Interest Rate')
ylim([0.01, 0.08])
title('Interest Rate Paths under CIR model with fixed step size')
set(gca, 'FontSize', 12 )
grid on; hold off;


figure; hold on;
for k = idx
    plot(tPaths{k}, PathsA{k}, 'LineWidth', 1);
end
xlim([0 T]);
xlabel('Time'); 
ylabel('Interest Rate');
title('Interest Rate Paths under CIR model with variable step size');
grid on; hold off;

%%
tp = tPaths{1};  rp = PathsA{1};
dt_seq = diff(tp);           % 每步的步长
t_mid  = tp(1:end-1);        % 步长对应的时刻（步的起点）

figure;
subplot(2,1,1);
plot(tp, rp, '-', 'LineWidth', 1.2); hold on;
yline(theta, 'k--', 'LineWidth', 1);
ylabel('Interest rate');
legend('r(t)', '\theta', 'Location','best');
title('Sample path and adaptive step size'); 
set(gca, 'FontSize', 12 );
grid on;

subplot(2,1,2);
stairs(t_mid, dt_seq, 'LineWidth', 1.2);
yline(hmax, 'r', 'LineWidth', 1);
xlabel('Time'); ylabel('Step size');
legend('Adaptive step size', 'h_{max}', 'Location','best'); 
set(gca, 'FontSize', 12 );
grid on;

%% 平均步数对比
tol_scan = logspace(-6, -2, 15);   % τ 从 1e-6 到 1e-2
navg = zeros(size(tol_scan));
for i = 1:numel(tol_scan)
    rng(seed);
    [~,~,tP] = CIR_Adap(r0,lambda,theta,eta,tol_scan(i),hmax,T,1e4,true);
    navg(i) = mean(cellfun(@numel,tP)) - 1;
end
figure;
semilogx(tol_scan, navg, '-o','LineWidth',1); hold on;
yline(T/hmax, 'k--', 'Fixed Step','LineWidth',1);
xlabel('Error Tolerance'); ylabel('Average Number of Time Steps');
title('Average Number of Time Steps versus Error Tolerance'); 
set(gca, 'FontSize', 12 );
grid on;

%%
alpha_near = 0.5;      % 近零区阈值:r < alpha_near*theta
alpha_hit  = 0.1;      % 触底阈值:  r < alpha_hit *theta

etas   = [0.05 0.07 0.10 0.14 0.18 0.22 0.26 0.283];
feller = 2*lambda*theta ./ etas.^2;
 
P_hit      = zeros(size(etas));   % P(min_t r < alpha_hit*theta)
frac_near0 = zeros(size(etas));   % 平均处于近零区的时间比例
 
N  = round(T/hmax);
dt = T/N;
 
for k = 1:numel(etas)
    eta = etas(k);
    rng(seed);                    % 每个 eta 用同一 seed,便于复现
    r    = r0*ones(M,1);
    rmin = r0*ones(M,1);
    near = zeros(M,1);
    for i = 1:N
        r_pos = max(r,0);                          % full truncation
        dW    = sqrt(dt)*randn(M,1);
        r     = r + lambda*(theta - r_pos)*dt + eta*sqrt(r_pos).*dW;
        rmin  = min(rmin, r);
        near  = near + (r < alpha_near*theta);
    end
    P_hit(k)      = mean(rmin < alpha_hit*theta);
    frac_near0(k) = mean(near)/N;
end

% ---- 作图:双 y 轴 ----
figure('Color','w'); 
yyaxis left
plot(feller, 100*P_hit, '-o','LineWidth',1.4,'MarkerSize',6); 
ylabel(sprintf('P(\\,r hits < %.2g\\theta\\,)  [%%]', alpha_hit));
ylim([0 100]);
yyaxis right
plot(feller, 100*frac_near0, '--s','LineWidth',1.4,'MarkerSize',6);
ylabel(sprintf('Time fraction with r < %.2g\\theta  [%%]', alpha_near));

set(gca,'XDir','reverse');        % 左侧 Feller 大(宽松)-> 右侧接近 1(收紧)
xlabel('Feller ratio  2\lambda\theta/\eta^2');
title('Near-zero behaviour of the CIR rate versus the Feller ratio');
grid on;
xline(1,'k:','Feller = 1','LabelVerticalAlignment','bottom');

% ---- 控制台汇总 ----
fprintf('%6s %8s %12s %14s\n','eta','Feller','P(r<0.1th)%','near0 frac %');
for k = 1:numel(etas)
    fprintf('%6.3f %8.2f %12.2f %14.2f\n', etas(k), feller(k), 100*P_hit(k), 100*frac_near0(k));
end
