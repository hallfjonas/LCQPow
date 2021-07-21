%% Run an edges problem
% Clean up and load libraries
clear all; clc; close all;

%addpath(fullfile(pwd, ".."));
addpath(fullfile(pwd, "../../../build/lib"))

%
%   min   (x-xk)'*Q*(x-xk)
%   s.t.   ||x||^2 = 1
%

QPO_dense = {};
QPO_sparse = {};
OSQP = {};

% General parameters
params.printLevel = 2;
params.stationarityTolerance = 1e-3;

Nvals = 10:20:200;
for j = 1:length(Nvals)
    % Circle discretization
    N = Nvals(j);
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
    params.x0 = ones(nz,1);
    params.x0(1:2) = xk;

    % Dense qpOASES
    params.qpSolver = 0;
    [xOpt,~,stats] = LCQPow(Q, g, L, R, [], [], [], [], A, lbA, ubA, params);
    QPO_dense.time(j) = stats.elapsed_time;
    QPO_dense.exit_flag(j) = stats.exit_flag;

    % Sparse methods:
    Q = sparse(Q);
    L = sparse(L);
    R = sparse(R);
    A = sparse(A);

    % Sparse qpOASES
    params.qpSolver = 1;
    [xOpt,~,stats] = LCQPow(Q, g, L, R, [], [], [], [], A, lbA, ubA, params);
    QPO_sparse.time(j) = stats.elapsed_time;
    QPO_sparse.exit_flag(j) = stats.exit_flag;

    params.qpSolver = 2;
    [xOpt,~,stats] = LCQPow(Q, g, L, R, [], [], [], [], A, lbA, ubA, params);
    OSQP.time(j) = stats.elapsed_time;
    OSQP.exit_flag(j) = stats.exit_flag;
end

%% Plots
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(1); hold on; box on; grid on;
plot(Nvals, QPO_dense.time, "DisplayName", "QPO DENSE");
plot(Nvals, QPO_sparse.time, "DisplayName",  "QPO SPARSE");
plot(Nvals, OSQP.time, "DisplayName", "OSQP");
legend;
set(gca,'yscale','log');

figure(2); hold on; box on; grid on;
plot(Nvals, QPO_dense.exit_flag, "DisplayName", "QPO DENSE");
plot(Nvals, QPO_sparse.exit_flag, "DisplayName",  "QPO SPARSE");
plot(Nvals, OSQP.exit_flag, "DisplayName", "OSQP");
legend;

