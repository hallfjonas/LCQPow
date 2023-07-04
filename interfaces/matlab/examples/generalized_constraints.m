%% Run a simple warm up problem
Q = [2 0; 0 2];
g = [-2; -2];
L = [-1 0];
R = [0 1];
lb_L = -2;
lb_R = 0;
ub_L = inf;
ub_R = inf;

% Run solver
[x, y, stats] = LCQPow(Q, g, L, R, lb_L, ub_L, lb_R, ub_R);

fprintf("x = [%g, %g]\n", x(1), x(2));

if (length(y) == 4)
    fprintf("y = [%g, %g, %g, %g]\n", y(1), y(2), y(3), y(4));
elseif (length(y) == 2)
    fprintf("y = [%g, %g]\n", y(1), y(2));
end