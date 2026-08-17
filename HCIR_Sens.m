%% Parameter Setting
S0 = 100; K =  95; T = 1; 
r0 = 0.04; eta = 0.15; lambda = 0.5; theta = 0.04; 
v0 = 0.04; kappa = 0.3; vbar = 0.04; gamma = 0.15;
rho_xv = -0.7; rho_xr = 0.6;
% sqrt(2*lambda*theta)

N = 200; % Number of fixed steps
M = 100000; % Number of simulations
tol = 0.00001;
hmax = T/200;

seed = 20260804;

%% Sensitivity of eta
eta_test = [0, 0.05, 0.1, 0.15, 0.2];
h_bench = 0.0005;
N_bench = round(T / h_bench);


% 0.2时刚好达到feller边界
[C_bench1, C_ref1, C_F1, C_V1, SE_bench, SE_F, SE_V, AvgSteps] = deal(zeros(1,length(eta_test)));


for i = 1:length(eta_test)
    if eta_test(i) == 0
        C_ref1(i) = Heston_Call(S0,K,r0,T,v0,kappa,vbar,gamma,rho_xv);
    else
        C_ref1(i) = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta_test(i), v0, kappa, vbar, gamma, rho_xv, rho_xr);
    end
    rng(seed);
    [ST_bench, I_bench, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta_test(i), v0, kappa, vbar, gamma, rho_xv, rho_xr, N_bench, M, T, false); 
    pb = exp(-I_bench) .* max(ST_bench-K, 0);
    C_bench1(i) = mean(pb); SE_bench(i) = std(pb)/sqrt(M); 

    rng(seed); 
    [ST_F, I_F, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta_test(i), v0, kappa, vbar, gamma, rho_xv, rho_xr, N, M, T, false); 
    pf = exp(-I_F) .* max(ST_F-K, 0);
    C_F1(i) = mean(pf); SE_F(i) = std(pf)/sqrt(M); 

    rng(seed); 
    [ST_V, I_V, ~, ~, ~, ns] = HCIR_Adaptive(r0, lambda, theta, eta_test(i), v0, kappa, vbar, gamma, tol, hmax, S0, rho_xv, rho_xr, M, T, false);
    pv = exp(-I_V) .* max(ST_V-K, 0);
    C_V1(i) = mean(pv); SE_V(i) = std(pv)/sqrt(M);

    fprintf('eta=%.3f: fixed steps = %.0f, avg steps = %.0f\n', eta_test(i), N, mean(ns));
    AvgSteps(i) = mean(ns);
end

Bias_F = C_F1 - C_bench1;
Bias_V = C_V1 - C_bench1;

AbsErr_F = abs(Bias_F);
AbsErr_V = abs(Bias_V);

SE_Bias_F = sqrt(SE_F.^2 + SE_bench.^2);
SE_Bias_V = sqrt(SE_V.^2 + SE_bench.^2);

figure; hold on;
errorbar(eta_test, C_F1, 1.96*SE_F, 'o-', 'LineWidth', 1);
plot(eta_test, C_ref1, 'o--', 'LineWidth', 1);
legend('Fixed-step Euler (95% CI)','Semi-analytical reference','Location', 'southeast');
xlabel('Interest-rate Volatility')
ylabel('European Call Option Price')
title('Effect of Interest-rate Volatility on Option Price under Fixed-step Scheme')
set(gca, 'FontSize', 12 )
grid on; hold off;

figure; hold on;
errorbar(eta_test, C_V1, 1.96*SE_V, 'o-', 'LineWidth', 1);
plot(eta_test, C_ref1, 'o--', 'LineWidth', 1);
legend('Adaptive Euler (95% CI)','Semi-analytical reference','Location', 'southeast');
xlabel('Interest-rate Volatility')
ylabel('European Call Option Price')
title('Effect of Interest-rate Volatility on Option Price under Adaptive-step Scheme')
set(gca, 'FontSize', 12 )
grid on; hold off;

%% Error vs eta
figure; hold on;
errorbar(eta_test, Bias_F, 1.96*SE_Bias_F, 'o-', 'LineWidth', 1.2);
errorbar(eta_test, Bias_V, 1.96*SE_Bias_V, 'o-', 'LineWidth', 1.2);
yline(0, 'k--', 'LineWidth', 1);

xlabel('Interest-rate volatility, \eta', 'Interpreter', 'tex');
ylabel('Price difference from fine-grid benchmark');
legend('Fixed-step Euler', 'Adaptive Euler', 'Zero error', ...
       'Location', 'best');
title('Discretisation Bias Relative to Fine-grid Benchmark');
grid on; hold off;

ApproxBias = C_ref1 - C_bench1;

figure; hold on;
plot(eta_test, ApproxBias, 'o-', 'LineWidth', 1);
xlabel('Interest-rate volatility');
ylabel('Price difference');
title('Price difference between reference and fine-grid benchmark');
grid on; hold off;

%%

figure; hold on;
plot(eta_test, AvgSteps, 'o-', 'LineWidth', 1);
yline(N, '--', 'Fixed-step Euler','LineWidth', 1.2);
xlabel('Interest-rate volatility, \eta', 'Interpreter', 'tex');
ylabel('Average number of adaptive time steps');
title('Effect of Interest-rate Volatility on Adaptive Time Steps');
legend('Adaptive Euler', 'Fixed-step Euler', 'Location', 'best');
set(gca, 'FontSize', 12);
grid on; hold off;


%%
nEta = length(eta_test);
R = 5;

Time_F = nan(R, nEta);
Time_V = nan(R, nEta);

for i = 1:nEta
    for j = 1:R

        rng(seed + j);

        tic;
        [ST_F, I_F, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta_test(i), ...
            v0, kappa, vbar, gamma, rho_xv, rho_xr, N, M, T, false);
        pf = exp(-I_F(:)) .* max(ST_F(:) - K, 0);
        C_temp_F = mean(pf);
        Time_F(j,i) = toc;

        rng(seed + j);

        tic;
        [ST_V, I_V, ~, ~, ~, ns] = HCIR_Adaptive( ...
            r0, lambda, theta, eta_test(i), v0, kappa, vbar, gamma, ...
            tol, hmax, S0, rho_xv, rho_xr, M, T, false);
        pv = exp(-I_V(:)) .* max(ST_V(:) - K, 0);
        C_temp_V = mean(pv); 
        Time_V(j,i) = toc;
    end
end

MeanTime_F = mean(Time_F, 1);
MeanTime_V = mean(Time_V, 1);

figure; hold on;
plot(eta_test, MeanTime_F, 'o-', 'LineWidth', 1.2);
plot(eta_test, MeanTime_V, 'o-', 'LineWidth', 1.2);

xlabel('Interest-rate volatility, \eta', 'Interpreter', 'tex');
ylabel('Average CPU time (seconds)');
title('Computational Cost as Interest-rate Volatility Changes');
legend('Fixed-step Euler', 'Adaptive Euler', 'Location', 'best');
grid on; hold off;

%% Sensitivity of rho_xr
eta_test2 = [0, 0.01, 0.05, 0.1];
rho_xr_max  = sqrt(1 - rho_xv^2);      
rho_xr_test = [0, 0.2, 0.4, 0.6, 0.7];
assert(all(rho_xr_test <= rho_xr_max), 'rho_xr exceeds positive-definite bound');

% When rho is large the log-price becomes the binding process, and the criterion cannot detect this.
nk = numel(eta_test2);  ni = numel(rho_xr_test);
C_ref_all = zeros(nk,ni);  C_F_all = zeros(nk,ni);  C_V_all = zeros(nk,ni);
SE_F_all  = zeros(nk,ni);  SE_V_all = zeros(nk,ni);  avg_steps = zeros(nk,ni);

for k = 1:nk
    eta = eta_test2(k);

    for i = 1:ni
        rxr = rho_xr_test(i);
        C_ref_all(k,i) = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rxr);
        
        rng(seed); 
        [ST_F, I_F, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rxr, N, M, T, false); 
        pf = exp(-I_F) .* max(ST_F-K, 0);
        C_F_all(k,i) = mean(pf); SE_F(i) = std(pf)/sqrt(M);  

        rng(seed); 
        [ST_V, I_V, ~, ~, ~, ns] = HCIR_Adaptive(r0, lambda, theta, eta, v0, kappa, vbar, gamma, tol, hmax, S0, rho_xv, rxr, M, T, false);
        pv = exp(-I_V) .* max(ST_V-K, 0);
        C_V_all(k,i) = mean(pv); SE_V(i) = std(pv)/sqrt(M);
        avg_steps(k,i) = mean(ns);
    end

    figure;hold on;
    ciF = 1.96*SE_F_all(k,:);
    fill([rho_xr_test, fliplr(rho_xr_test)], ...
         [C_F_all(k,:)+ciF, fliplr(C_F_all(k,:)-ciF)], ...
         [0 0.45 0.74], 'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');

    plot(rho_xr_test, C_F_all(k,:),  '-o','LineWidth',1)
    plot(rho_xr_test, C_V_all(k,:),  '-o','LineWidth',1)
    plot(rho_xr_test, C_ref_all(k,:),'-','LineWidth',1.5)
    xlabel('Correlation Coefficient, \rho_{xr}', 'Interpreter','tex')
    ylabel('Option price'); ylim([12.9,13.15]);
    lgd = legend('Fixed-step Euler','Adaptive Euler','Reference price');
    title(sprintf('\\eta = %.2f', eta), 'Interpreter','tex')
    set(gca, 'FontSize', 12 )
    grid on; hold off;
end

%%
eta_check = [0.05, 0.10];
nEta = numel(eta_check);
nRho = numel(rho_xr_test);

% Fine-grid setting: ten times finer than the ordinary fixed-step grid
h_bench = 0.0005;
N_bench = round(T / h_bench);

[C_bench2, C_F2, C_V2, SE_bench2, SE_F2, SE_V2, AvgSteps2] = deal(zeros(nEta, nRho));

for j = 1:nEta
    eta_now = eta_check(j);
    for i = 1:nRho
        rho_now = rho_xr_test(i);

        % Fine-grid benchmark: use an independent seed
        rng(seed + j + i);
        [ST_b, I_b, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta_now, ...
            v0, kappa, vbar, gamma, rho_xv, rho_now, ...
            N_bench, M, T, false);
        pb = exp(-I_b(:)) .* max(ST_b(:) - K, 0);
        C_bench2(j,i) = mean(pb);
        SE_bench2(j,i) = std(pb) / sqrt(numel(pb));

         % Ordinary fixed-step scheme: independent seed
        rng(seed + j + i);
        [ST_f, I_f, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta_now, ...
            v0, kappa, vbar, gamma, rho_xv, rho_now, ...
            N, M, T, false);
        pf = exp(-I_f(:)) .* max(ST_f(:) - K, 0);
        C_F2(j,i) = mean(pf);
        SE_F2(j,i) = std(pf) / sqrt(numel(pf));

        % Adaptive-step scheme: independent seed
        rng(seed +j + i);
        [ST_a, I_a, ~, ~, ~, ~] = HCIR_Adaptive( ...
            r0, lambda, theta, eta_now, ...
            v0, kappa, vbar, gamma, tol, hmax, ...
            S0, rho_xv, rho_now, M, T, false);

        pa = exp(-I_a(:)) .* max(ST_a(:) - K, 0);
        C_V2(j,i) = mean(pa);
        SE_V2(j,i) = std(pa) / sqrt(numel(pa));
    end
end

% Signed price differences from the fine-grid benchmark
Diff_fixed = C_F2 - C_bench2;
Diff_adapt = C_V2 - C_bench2;

SEdiff_fixed = sqrt(SE_F2.^2 + SE_bench2.^2);
SEdiff_adapt = sqrt(SE_V2.^2 + SE_bench2.^2);

% Plot: 95% confidence intervals for price differences
figure;

for j = 1:nEta
    figure; hold on;

    errorbar(rho_xr_test, Diff_fixed(j,:), 1.96*SEdiff_fixed(j,:), ...
        'o-', 'LineWidth', 1, 'MarkerSize', 6);
    errorbar(rho_xr_test, Diff_adapt(j,:), 1.96*SEdiff_adapt(j,:), ...
        'o-', 'LineWidth', 1, 'MarkerSize', 6);
    yline(0, 'k--', 'LineWidth', 1);

    xlabel('Correlation coefficient');
    ylabel('Price difference from fine-grid benchmark');
    title(sprintf('\\eta = %.2f', eta_check(j)));
    legend('Fixed-step Euler', 'Adaptive Euler', 'Zero difference', ...
        'Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);
    hold off;
end





%% Sensitivity of T
T_test = [1, 2, 5, 10, 20];
[C_ref, C_F, C_V,C_bench_T,SE_bench_T, SE_F, SE_V, nsF, nsV] = deal(zeros(1,length(T_test)));
h_fixed = 0.05;

for i = 1:length(T_test)
    rng(seed+i);
    [ST_bench, I_bench, ~] = HCIR_Fixed( ...
        S0, r0, lambda, theta, eta, ...
        v0, kappa, vbar, gamma, rho_xv, rho_xr, ...
        T_test(i)/h_fixed*10, M, T_test(i), false);

    pb = exp(-I_bench(:)) .* max(ST_bench(:) - K, 0);
    C_bench_T(i) = mean(pb);
    SE_bench_T(i) = std(pb) / sqrt(numel(pb));

    rng(seed+i); 
    [ST_F, I_F, ~] = HCIR_Fixed(S0, r0, lambda, theta, eta, v0, kappa, vbar, gamma, rho_xv, rho_xr, T_test(i)/h_fixed, M, T_test(i), false); 
    pf = exp(-I_F) .* max(ST_F-K, 0);
    C_F(i) = mean(pf); SE_F(i) = std(pf)/sqrt(M); 

    rng(seed+i); 
    [ST_V, I_V, ~, ~, ~, ns] = HCIR_Adaptive(r0, lambda, theta, eta, v0, kappa, vbar, gamma, tol, hmax, S0, rho_xv, rho_xr, M, T_test(i), false);
    pv = exp(-I_V) .* max(ST_V-K, 0);
    C_V(i) = mean(pv); SE_V(i) = std(pv)/sqrt(M);
end

Diff_F_T = C_F - C_bench_T;
Diff_V_T = C_V - C_bench_T;
SE_Diff_F_T = sqrt(SE_F.^2 + SE_bench_T.^2);
SE_Diff_V_T = sqrt(SE_V.^2 + SE_bench_T.^2);

figure; hold on;
plot(T_test,C_F,'LineWidth',1)
plot(T_test,C_V,'LineWidth',1)
plot(T_test,C_bench_T,'LineWidth',1)
legend('Fixed-step Euler', 'Adaptive Euler', 'fine grid Price');
xlabel('Time to Maturity, T')
ylabel('Option price')
title('Effect of Time to Maturity on European Call Option Price')
set(gca, 'FontSize', 12 )
grid on; hold off;

% option price随着maturity变长而上升，期限越长时间价值越大

figure; hold on;

errorbar(T_test, Diff_F_T, 1.96 * SE_Diff_F_T, ...
    'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

errorbar(T_test, Diff_V_T, 1.96 * SE_Diff_V_T, ...
    'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

yline(0, 'k--', 'LineWidth', 1.1);

xlabel('Time to maturity, T');
ylabel('Price difference from fine-grid benchmark');
title('Pricing Differences Relative to Fine-grid Benchmark');
legend('Fixed-step Euler', 'Adaptive Euler', 'Zero difference', ...
    'Location', 'best');

set(gca, 'FontSize', 12);
grid on;
hold off;


%% Cost comparison: fixed h versus adaptive tolerance
T_cost = 5;

h_test = [0.05, 0.02, 0.01, 0.005];
delta_test = [0.01, 0.005, 0.001, 0.0005, 0.0001];

R = 1;          % repeated timing runs
M_cost = M;     % same number of Monte Carlo paths

nH = numel(h_test);
nD = numel(delta_test);

% One row for each (h, delta) combination
nRows = nH * nD;

H_col = nan(nRows,1);
Delta_col = nan(nRows,1);
FixedSteps_col = nan(nRows,1);
FixedCPU_col = nan(nRows,1);
AdaptiveSteps_col = nan(nRows,1);
AdaptiveCPU_col = nan(nRows,1);

row = 0;

for ih = 1:nH

    h_now = h_test(ih);
    N_now = round(T_cost / h_now);

    % Fixed-step CPU time for this h
    FixedRuns = nan(R,1);

    for j = 1:R
        rng(seed + 10000*ih + j);

        tic;
        [ST_F, I_F, ~] = HCIR_Fixed( ...
            S0, r0, lambda, theta, eta, ...
            v0, kappa, vbar, gamma, rho_xv, rho_xr, ...
            N_now, M_cost, T_cost, false);

        payoffF = exp(-I_F(:)) .* max(ST_F(:) - K, 0);
        mean(payoffF); %#ok<VUNUS>
        FixedRuns(j) = toc;
    end

    FixedCPU_now = median(FixedRuns);

    % Adaptive-step CPU time for each tolerance
    for id = 1:nD

        delta_now = delta_test(id);

        AdaptiveRuns = nan(R,1);
        AdaptiveStepRuns = nan(R,1);

        for j = 1:R
            rng(seed + 20000*ih + 1000*id + j);

            tic;
            [ST_V, I_V, ~, ~, ~, ns] = HCIR_Adaptive( ...
                r0, lambda, theta, eta, ...
                v0, kappa, vbar, gamma, ...
                delta_now, h_now, S0, ...
                rho_xv, rho_xr, M_cost, T_cost, false);

            payoffV = exp(-I_V(:)) .* max(ST_V(:) - K, 0);
            mean(payoffV); %#ok<VUNUS>

            AdaptiveRuns(j) = toc;
            AdaptiveStepRuns(j) = mean(ns);
        end

        row = row + 1;

        H_col(row) = h_now;
        Delta_col(row) = delta_now;
        FixedSteps_col(row) = N_now;
        FixedCPU_col(row) = FixedCPU_now;
        AdaptiveSteps_col(row) = mean(AdaptiveStepRuns);
        AdaptiveCPU_col(row) = median(AdaptiveRuns);
    end
end

% Output table
CostTable = table(H_col, Delta_col, FixedSteps_col, FixedCPU_col, ...
    AdaptiveSteps_col, AdaptiveCPU_col, ...
    'VariableNames', {'h','delta','FixedSteps','FixedCPUTime', ...
                      'AdaptiveSteps','AdaptiveCPUTime'});

disp(CostTable);


%% Computational-cost analysis
% Assumes all model parameters already exist in the workspace.

R = 1;                 % Number of repeated timing runs
M_cost = M;            % Keep the same number of paths as pricing experiments
h_fixed = 0.05;
T_base = 5;
tol = 0.0001;

P = struct( ...
    'S0', S0, 'K', K, 'r0', r0, ...
    'lambda', lambda, 'theta', theta, 'eta', eta, ...
    'v0', v0, 'kappa', kappa, 'vbar', vbar, 'gamma', gamma, ...
    'rho_xv', rho_xv, 'rho_xr', rho_xr, ...
    'M', M_cost, 'seed', seed);

% Panel A: Different fixed time steps, T = 1
N_test = [50, 100, 200, 500, 1000];

TimeF_A = nan(size(N_test));
TimeV_A = nan(size(N_test));
StepsV_A = nan(size(N_test));

for i = 1:numel(N_test)
    N_now = N_test(i);
    hmax_now = T_base / N_now;  % Compare with the same maximum step

    [TimeF_A(i), TimeV_A(i), StepsV_A(i)] = ...
        measure_cost(P, T_base, N_now, hmax_now, tol, R);
end

PanelA = table(N_test(:), N_test(:), TimeF_A(:), StepsV_A(:), TimeV_A(:), ...
    'VariableNames', {'Setting_N','FixedSteps','FixedCPUTime', ...
                      'AdaptiveSteps','AdaptiveCPUTime'});

disp('Panel A: Different numbers of fixed time steps (T = 5)');
disp(PanelA);

%% Panel B: Different maturities, with h_fixed held constant
T_test = [1, 2, 5, 10, 20];

TimeF_B = nan(size(T_test));
TimeV_B = nan(size(T_test));
StepsF_B = nan(size(T_test));
StepsV_B = nan(size(T_test));

for i = 1:numel(T_test)
    T_now = T_test(i);
    N_now = round(T_now / h_fixed);

    StepsF_B(i) = N_now;

    [TimeF_B(i), TimeV_B(i), StepsV_B(i)] = ...
        measure_cost(P, T_now, N_now, h_fixed, tol, R);
end

PanelB = table(T_test(:), StepsF_B(:), TimeF_B(:), StepsV_B(:), TimeV_B(:), ...
    'VariableNames', {'Setting_T','FixedSteps','FixedCPUTime', ...
                      'AdaptiveSteps','AdaptiveCPUTime'});

disp('Panel B: Different maturities (h_fixed = 0.05)');
disp(PanelB);

%% Panel C: Different adaptive tolerances, T = 5
tol_test = [0.01, 0.005, 0.001, 0.0005, 0.0001];

N_base = round(T_base / h_fixed);

TimeF_C = nan(size(tol_test));
TimeV_C = nan(size(tol_test));
StepsV_C = nan(size(tol_test));

for i = 1:numel(tol_test)
    [TimeF_C(i), TimeV_C(i), StepsV_C(i)] = ...
        measure_cost(P, T_base, N_base, h_fixed, tol_test(i), R);
end

PanelC = table(tol_test(:), repmat(N_base,numel(tol_test),1), ...
    TimeF_C(:), StepsV_C(:), TimeV_C(:), ...
    'VariableNames', {'SettingTolerance','FixedSteps','FixedCPUTime', ...
                      'AdaptiveSteps','AdaptiveCPUTime'});

% Display results

disp('Panel C: Different tolerances (T = 5, h_fixed = 0.05)');
disp(PanelC);