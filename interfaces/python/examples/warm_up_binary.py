import lcqpow
import numpy as np


print("Preparing dense warm up problem...")

# Setup data of first QP.
H = np.array([[2.0, 0.0], [0.0, 2.0]])
g = np.array([-2.0, -2.0])

# 0 <= x _|_ y >= 0
# 0 <= x _|_ 0.5 - x >= 0
S1 = np.array([[1.0, 1.0], [0.0, 0.0]])
S2 = np.array([[0.0, -1.0], [1.0, 0.0]])
lbS1 = np.array([0.0, 0.0])
lbS2 = np.array([0.0, -0.5])

x0 = np.array([0.0, 0.0])

nV = 2 
nC = 0
nComp = 2
lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)

options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.QPOASES_DENSE)
lcqp.setOptions(options)

# Solve first LCQP
retVal = lcqp.loadLCQP(H=H, g=g, S1=S1, S2=S2, lbS1=lbS1, lbS2=lbS2, x0=x0)
if retVal is not lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal is not lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to solve LCQP.")


# Solve another LCQP
x0[0] = 0.0
x0[1] = 3000.0
options.setSolveZeroPenaltyFirst(False)
options.setInitialPenaltyParameter(10.0)
lcqp.setOptions(options)
retVal = lcqp.loadLCQP(H=H, g=g, S1=S1, S2=S2, lbS1=lbS1, lbS2=lbS2, x0=x0)
if retVal is not lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal is not lcqpow.ReturnValue.SUCCESSFUL_RETURN:
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