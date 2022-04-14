import lcqpow
import numpy as np


print("Preparing unit circle optimization problem...")

# Set dimensions
N = 100
nV = 2 + 2*N
nC = N+1
nComp = N

# Reference 
x_ref= np.array([0.5, -0.6])

# Set up LCQP object
lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)
options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.OSQP_SPARSE)
options.setStationarityTolerance(10e-3)
lcqp.setOptions(options)

# Allocate vectors
Q = np.zeros([nV, nV])
g = np.zeros(nV)
L = np.zeros([nComp, nV])
R = np.zeros([nComp, nV])
A = np.zeros([nC, nV])
lbA = np.zeros(nC)
ubA = np.zeros(nC)
x0 = np.zeros(nV)

x0[0] = x_ref[0]
x0[1] = x_ref[1]

# Assign problem data
Q[0, 0] = 17 
Q[1, 1] = 17
Q[0, 1] = -15
Q[1, 0] = -15

# Regularization on Q
for i in range(2, nV, 1):
    Q[i, i] = 5e-12

# Objective linear term
Qx = np.array([[17, -15], [-15, 17]])
g[:2] = - np.dot(Qx, x_ref)

# Constraints
for i in range(N):
    # Equality constraint ([cos sin]'*x + lambda = 1)
    A[i, 0] = np.cos((2*np.pi*i)/N)
    A[i, 1] = np.sin((2*np.pi*i)/N)
    A[i, 2+2*i] = 1.0

    # Convex combination constraint (sum theta = 1)
    A[N, 3+2*i] = 1.0

    # Complementarity constraints (lamba*theta = 0)
    L[i, 2+2*i] = 1.0
    R[i, 3+2*i] = 1.0

    # constraint bounds
    lbA[i] = 1.0
    ubA[i] = 1.0

    x0[2*i+2] = 1.0
    x0[2*i+3] = 1.0

# Constraint bound for last constraint
lbA[N] = 1.0
ubA[N] = 1.0

# Solve first LCQP
retVal = lcqp.loadLCQP(Q=Q, g=g, L=L.T, R=R.T, A=A.T, lbA=lbA, ubA=ubA, x0=x0)
if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to solve LCQP.")

stats = lcqpow.OutputStatistics()
xOpt = lcqp.getPrimalSolution()
yOpt = lcqp.getDualSolution()
lcqp.getOutputStatistics(stats)
print("xOpt = ", xOpt[0:2])
print("i = ", stats.getIterTotal())
print("k = ", stats.getIterOuter())
print("rho = ", stats.getRhoOpt())
print("WSR = ", stats.getSubproblemIter())