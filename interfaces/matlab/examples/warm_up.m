%% Run a simple warm up problem
% Clean up and load libraries
clear all; clc; close all;

% Path to mex function (c++ solver)
addpath(fullfile(pwd, ".."));

% Simple QP
Q = [2 0; 0 2];
g = [-2; -2];
L = [1 0];
R = [0 1];
lb = [-inf; -inf];
ub = [inf; inf];

% Algorithm parameters (C++)
params.x0 = [1; 1];
params.initialPenaltyParameter = 0.01;
params.penaltyUpdateFactor = 2;

% Run solver
[x, y] = LCQPanther(Q, g, L, R, lb, ub, params);

fprintf("x = [%g, %g]\n", x(1), x(2));
fprintf("y = [%g, %g, %g, %g]\n", y(1), y(2), y(3), y(4));
