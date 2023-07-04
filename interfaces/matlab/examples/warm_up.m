%% Run a simple warm up problem
Q = [2 0; 0 2];
g = [-2; -2];
L = [1 0];
R = [0 1];

% Algorithm parameters (C++)
params.initialPenaltyParameter = 0.01;
params.penaltyUpdateFactor = 2;
params.qpSolver = 1;

% Pass qpOASES options/ OSQP options like this:
%    REMARK: qpOASES and osqp-matlab interfaces must be accessible
% params.qpOASES_options = qpOASES_options();
% params.OSQP_options = osqp.default_settings();

if params.qpSolver > 0
    Q = sparse(Q);
    L = sparse(L);
    R = sparse(R);
end

% Run solver
[x, y] = LCQPow(Q, g, L, R, [], [], [], [], params);

fprintf("x = [%g, %g]\n", x(1), x(2));

if (length(y) == 4)
    fprintf("y = [%g, %g, %g, %g]\n", y(1), y(2), y(3), y(4));
elseif (length(y) == 2)
    fprintf("y = [%g, %g]\n", y(1), y(2));
end
