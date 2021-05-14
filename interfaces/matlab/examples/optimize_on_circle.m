%% Run an edges problem
% Clean up and load libraries
clear all; clc; close all;

addpath(fullfile(pwd, ".."));

%
%   min   (x-xk)'*Q*(x-xk)
%   s.t.   ||x||^2 = 1
%

% Circle discretization
N = 5;
nz = 2 + 2*N;

Qx = [17 -15; -15 17];
xk = [0.5; -0.6];
g = zeros(nz,1);

Q = eps*eye(nz);
Q(1:2,1:2) = Qx;
g(1:2) = -xk'*Qx;

obj = @(x) (x(1:2) - xk)'*Qx*(x(1:2) - xk);

A = zeros(N+1, nz);
lbA = ones(N+1, 1);
ubA = lbA;

L = zeros(N, nz);
R = zeros(N, nz);

lb = zeros(nz,1); lb(1:2) = -inf(2,1);
ub = inf(nz,1);

for i = 1:N
    % Equality constraint ([cos sin]'*x + lambda = 1)
    A(i, 1:2) = [cos(2*pi*i/N) sin(2*pi*i/N)];
    A(i, 2 + 2*i - 1) = 1;

    % Convex combination constraint (sum theta = 1)
    A(N+1, 2+2*i) = 1;

    % Complementarity constraint
    L(i, 2 + 2*i - 1) = 1;
    R(i, 2 + 2*i) = 1;
end

% Algorithm parameters
params.x0 = ones(nz,1); params.x0(1:2) = xk;
params.initialPenaltyParameter = 0.001;
params.penaltyUpdateFactor = 2;
params.solveZeroPenaltyFirst = false;
params.printLevel = 2;

[xOpt,~,stats] = LCQPanther(Q, g, L, R, A, lbA, ubA, lb, ub, params);

% Set to latex (only relevant for plotting)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Performance test
% Call solver for various rho init values
rho0 = [1.5*eps];
while(rho0(end) < 10^2)
    rho0 = [rho0, rho0(end)*1.5];
end

time_vals = [];
obj_vals = [];
iter_vals = [];
outer_vals = [];
for i = 1:length(rho0)
    params.initialPenaltyParameter = rho0(i);

    tic;
    [xOpt,~,stats] = LCQPanther(Q, g, L, R, A, lbA, ubA, lb, ub, params);

    % Evaluation
    time_vals = [time_vals, toc];
    obj_vals = [obj_vals, obj(xOpt(:,end))];
    iter_vals = [iter_vals, stats.iters_total];
    outer_vals = [outer_vals, stats.iters_outer];
end

%% Performance Plot
CreatePerformancePlot(rho0, time_vals, obj_vals, iter_vals, outer_vals);

%% Step Plot
params.storeSteps = true;
params.printStats = true;

rho0 = [1, 10, 12.4094];
CreateEdgesPlot(xk, N, Qx);

% global minimum
params.R0 = rho0(1);
xOpt_1 = LCQPanther(Q, g, A, lb, ub, lbA, ubA, L, R, params);
AddSteps(xOpt_1, "$\rho_0 = 1$", 'b', 'x');

% (strict) local minimum
params.R0 = rho0(2);
xOpt_10 = LCQPanther(Q, g, A, lb, ub, lbA, ubA, L, R, params);
AddSteps(xOpt_10, "$\rho_0 = 10$", 'r', 'x');

% Uncomment to plot this M-stationary point
%params.R0 = rho0(3);
%xOpt_12 = LCQP(Q, g, A, lb, ub, lbA, ubA, L, R, params);
%AddSteps(xOpt_12, "$\rho_0 = 12.4094$", 'g', 'o');

saveas(gcf,strcat('steps', '.png'));

%% Objective versus angle
ObjectiveVersusAngle(obj);