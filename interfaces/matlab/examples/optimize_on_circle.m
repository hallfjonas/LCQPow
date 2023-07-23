%% Build the LCQP
%
%   min   (x-xk)'*Q*(x-xk)
%   s.t.   ||x||^2 = 1
%

% Circle discretization
N = 200;
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
params.x0 = ones(nz,1); params.x0(1:2) = xk;
params.initialPenaltyParameter = 0.001;
params.penaltyUpdateFactor = 2;
params.solveZeroPenaltyFirst = true;
params.printLevel = 0;

% Use sparse solver
params.qpSolver = 1;
if (~issparse(Q))
    Q = sparse(Q);
    L = sparse(L);
    R = sparse(R);
    A = sparse(A);
end

% Set to latex (only relevant for plotting)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Performance test
% Call solver for various rho init values
rho0 = [10e-8];
while(rho0(end) < 10^1)
    rho0 = [rho0, rho0(end)*1.5];
end

c_solver = struct;
m_solver = struct;

plot_maps = containers.Map;

c_solver.time_vals = zeros(length(rho0),1);
c_solver.obj_vals = zeros(length(rho0),1);
c_solver.iter_vals = zeros(length(rho0),1);
c_solver.outer_vals = zeros(length(rho0),1);

for i = 1:length(rho0)
    params.initialPenaltyParameter = rho0(i);

    % Run c solver
    tic;
    [xOpt,~,stats] = LCQPow(Q, g, L, R, [], [], [], [], A, lbA, ubA, params);

    % Evaluation
    c_solver.time_vals(i) = toc;
    c_solver.obj_vals(i) = obj(xOpt(:,end));
    c_solver.iter_vals(i) = stats.iters_total;
    c_solver.outer_vals(i) = stats.iters_outer;
end

plot_maps("C++") = c_solver;

% Performance Plot
CreatePerformancePlot(rho0, plot_maps);

% Objective versus angle
ObjectiveVersusAngle(obj);

%% Plot Functions
function[] = CreatePerformancePlot(rho, plot_maps)

names = plot_maps.keys;

figure;
subplot(2,2,1); grid on; box on;
for i = 1:length(names)
    plot_struct = plot_maps(names{i});
    loglog(rho, plot_struct.time_vals); hold on;
end
ylabel("CPU time in [s]");
xlabel("Initial penalty value");
legend(names);

subplot(2,2,2); grid on; box on;
for i = 1:length(names)
    plot_struct = plot_maps(names{i});
    loglog(rho, plot_struct.obj_vals); hold on;
end
ylabel("Objective values");
xlabel("Initial penalty value");
legend(names);

subplot(2,2,3); grid on; box on;
for i = 1:length(names)
    plot_struct = plot_maps(names{i});
    loglog(rho, plot_struct.iter_vals); hold on;
end
ylabel("Number of iterations");
xlabel("Initial penalty value");
legend(names);

subplot(2,2,4); grid on; box on;
for i = 1:length(names)
    plot_struct = plot_maps(names{i});
    loglog(rho, plot_struct.outer_vals); hold on;
end
ylabel("Number of outer iterations");
xlabel("Initial penalty value");
legend(names);

end

% Objective versus angle
function[] = ObjectiveVersusAngle(obj)

figure; hold on; box on;

% Number of linear segments approximating the circle
N = 100;

obj_vals = zeros(N, 1);
angle_vals = zeros(N, 1);

% Fill objective values
vertices = [cos(2*pi*(1:N)/N); sin(2*pi*(1:N)/N)];
for i = 1:N
    angle_vals(i) = 2*pi*i/N;
    obj_vals(i) = obj(vertices(:,i));
end

% Create the plot
plot(angle_vals, obj_vals + 10);
xlabel("angle of $x$");
ylabel("objective value");

% Set to log scale
set(gca,'yscale','log')

end
