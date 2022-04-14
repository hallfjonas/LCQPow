import lcqpow
import numpy as np


print("Preparing OSQP warm up problem...")

# Setup data of first QP.
Q_data = np.array([2.0, 2.0])
Q_i = [0, 1]
Q_p = [0, 1, 2]

g = np.array([-2.0, -2.0])

L_data = np.array([1.0])
L_i = [0]
L_p = [0, 1, 1]

R_data = np.array([1.0])
R_i = [0]
R_p = [0, 0, 1]

nV = 2 
nC = 0
nComp = 1

Q = lcqpow.cscWrapper(nV, nV, Q_p[nV], Q_data, Q_i, Q_p)
L = lcqpow.cscWrapper(nComp, nV, L_p[nV], L_data, L_i, L_p)
R = lcqpow.cscWrapper(nComp, nV, R_p[nV], R_data, R_i, R_p)

lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)

options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.OSQP_SPARSE)
lcqp.setOptions(options)

retVal = lcqp.loadLCQP(Q=Q, g=g, L=L, R=R)
if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to solve LCQP.")

stats = lcqpow.OutputStatistics()
xOpt = lcqp.getPrimalSolution()
yOpt = lcqp.getDualSolution()
lcqp.getOutputStatistics(stats)
print("xOpt = ", xOpt)
print("yOpt = ", yOpt)
print("i = ", stats.getIterTotal())
print("k = ", stats.getIterOuter())
print("rho = ", stats.getRhoOpt())
print("WSR = ", stats.getSubproblemIter())