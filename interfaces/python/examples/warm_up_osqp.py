import lcqpow
import numpy as np


print("Preparing OSQP warm up problem...")

# Setup data of first QP.
H_data = np.array([2.0, 2.0])
H_i = [0, 1]
H_p = [0, 1, 2]

g = np.array([-2.0, -2.0])

S1_data = np.array([1.0])
S1_i = [0]
S1_p = [0, 1, 1]

S2_data = np.array([1.0])
S2_i = [0]
S2_p = [0, 0, 1]

nV = 2 
nC = 0
nComp = 1

H = lcqpow.cscWrapper(nV, nV, H_p[nV], H_data, H_i, H_p)
S1 = lcqpow.cscWrapper(nComp, nV, S1_p[nV], S1_data, S1_i, S1_p)
S2 = lcqpow.cscWrapper(nComp, nV, S2_p[nV], S2_data, S2_i, S2_p)

lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)

options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.OSQP_SPARSE)
lcqp.setOptions(options)

retVal = lcqp.loadLCQP(H=H, g=g, S1=S1, S2=S2)
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