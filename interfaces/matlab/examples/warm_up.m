%% Run a simple warm up problem
% Clean up and load libraries
clear all; clc; close all;

% Path to mex function (c++ solver)
addpath(fullfile(pwd, ".."));

% Path to LCQP and qpOASES (matlab solver)
addpath("~/LCQP");
addpath("~/qpOASES/interfaces/matlab");

% Simple QP
Q = [2 0; 0 2];
g = [-2; -2];
L = [1 0];
R = [0 1];
A = [1 -1];
lbA = [-0.5];
ubA = [inf];
lb = [-inf; -inf];
ub = [inf; inf];

% Algorithm parameters (C++)
params.x0 = [1; 1];
params.initialPenaltyParameter = 0.01;
params.penaltyUpdateFactor = 2;
params.solveZeroPenaltyFirst = true;
params.printLevel = 0;

% Run solver
nexp = 100;
times_c = zeros(nexp, 1);
times_m = zeros(nexp, 1);
for i = 1:nexp
    tic;
    [xopt, yopt, stats] = LCQPanther(Q, g, L, R, A, lbA, ubA, lb, ub, params);
    times_c(i) = toc;

    tic;
    LCQP(Q, g, L, R, A, lbA, ubA, lb, ub, params);
    times_m(i) = toc;
end

% Plot results
figure;
subplot(1,2,1); hold on; box on;
semilogy(times_c, 'r--', 'DisplayName', 'C++');
semilogy(times_m, 'b--', 'DisplayName', 'MATLAB');

subplot(1,2,2); hold on; box on;
plot(times_m./times_c, '-', 'DisplayName', 'Speed Factor');
plot(1:nexp, repmat(mean(times_m./times_c), nexp, 1), '-', 'DisplayName', 'Speed Factor');
