% Simple QP
Q = [2*100^2 0; 0 2];
g = [-200; -2];
L = [1 0];
R = [0 1];

% Algorithm parameters (C++)
params.x0 = [1; 1];
params.initialPenaltyParameter = 10;
params.penaltyUpdateFactor = 2;
params.solveZeroPenaltyFirst = true;
params.printLevel = 2;

% Run solver
nexp = 1;
for i = 1:nexp
    tic;
    [xopt, yopt, stats] = LCQPow(Q, g, L, R, [], [], [], [], params);
    times_c(i) = toc;

    %tic;
    %LCQP(Q, g, A, lb, ub, lbA, ubA, L, R, params_mat)
    %times_m(i) = toc;
end

