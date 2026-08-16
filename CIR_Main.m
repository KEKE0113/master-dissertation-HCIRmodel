% Parameters
S0 = 100;
K =  95;

r0 = 0.04;
lambda = 0.5;
theta = 0.04; 
eta = 0.15;

N = 200;    % Number of fixed steps
M = 100000; % Number of simulations
T = 1; 

tol = 0.00001;
hmax = T/200;

seed = 20260731;

% CIR simulation
rng(seed);  [rF, PathsF] = CIR_Fixed(r0, lambda, theta, eta, T, N, M, true);
rng(seed);  [rA, PathsA, tPaths, nsteps] = CIR_Adap(r0, lambda, theta, eta, tol, hmax, T, M, true);

% Tolerance experiment

tol_scan = logspace(-6, -2, 15);
navg = zeros(size(tol_scan));

for i = 1:numel(tol_scan)
    rng(seed);

    [~, ~, ~, nsteps_tol] = CIR_Adap(r0, lambda, theta, eta, tol_scan(i), hmax, T, 10000, false);

    navg(i) = mean(nsteps_tol);
end
