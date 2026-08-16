%% Fixed-step CIR paths (Figure 5)

t = linspace(0, T, N + 1);

num_paths = 20;
rng(seed);
idx = randperm(M, num_paths);

figure;
plot(t, PathsF(idx, :)', 'LineWidth', 1);

xlabel('Time');
ylabel('Short Rate');
title('CIR Rate Paths: Fixed-Step Scheme');

set(gca, 'FontSize', 12);
grid on;

%% Adaptive CIR paths (Figure 5)

figure;
hold on;

for k = idx
    plot(tPaths{k}, PathsA{k}, 'LineWidth', 1);
end

xlim([0, T]);
xlabel('Time');
ylabel('Short Rate');
title('CIR Rate Paths: Adaptive Scheme');

set(gca, 'FontSize', 12);
grid on;
hold off;

%% Adaptive step-size behaviour (Figure 6)

tp = tPaths{1};  
rp = PathsA{1};

dt_seq = diff(tp);          
t_mid  = tp(1:end-1);      

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

xlabel('Time'); 
ylabel('Step size');
legend('Adaptive step size', 'h_{max}', 'Location','best'); 
set(gca, 'FontSize', 12 );
grid on;

%% Average number of adaptive steps (Figure 7)

tol_scan = logspace(-6, -2, 15);   % tolerance selected from 1e-6 to 1e-2
navg = zeros(size(tol_scan));

for i = 1:numel(tol_scan)
    rng(seed);

    [~, ~, ~, nsteps_tol] = CIR_Adap(r0, lambda, theta, eta, tol_scan(i), hmax, T, 1e4, false);

    navg(i) = mean(nsteps_tol);
end

figure;
semilogx(tol_scan, navg, '-o','LineWidth',1); hold on;

yline(T/hmax, 'k--', 'Fixed Step','LineWidth',1);

xlabel('Error Tolerance'); 
ylabel('Average Number of Time Steps');
title('Average Number of Time Steps versus Error Tolerance'); 

set(gca, 'FontSize', 12 );
grid on; hold off;
