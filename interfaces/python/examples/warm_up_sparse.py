import lcqpow
import numpy as np


print("Preparing sparse warm up problem...")

# Setup data of first QP.
Q = np.array([[2.0, 0.0], 
              [0.0, 2.0]])
g = np.array([-2.0, -2.0])
L = np.array([[1.0, 0.0]])
R = np.array([[0.0, 1.0]])

nV = 2 
nC = 0
nComp = 1
lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)

options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.QPOASES_SPARSE)
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