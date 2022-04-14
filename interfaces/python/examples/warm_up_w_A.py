import lcqpow
import numpy as np


print("Preparing dense warm up problem...")

# Setup data of first QP.
Q = np.array([[2.0, 0.0], 
              [0.0, 2.0]])
g = np.array([-2.0, -2.0])
L = np.array([[1.0, 0.0]])
R = np.array([[0.0, 1.0]])
A = np.array([[1.0, -1.0]])
lbA = np.array([-0.5])
ubA = np.array([np.infty])

nV = 2 
nC = 1
nComp = 1
lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)

options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
lcqp.setOptions(options)

retVal = lcqp.loadLCQP(Q=Q, g=g, L=L, R=R, A=A, lbA=lbA, ubA=ubA)
if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to solve LCQP.")

xOpt = lcqp.getPrimalSolution()
yOpt = lcqp.getDualSolution()
print("xOpt = ", xOpt)
print("yOpt = ", yOpt)