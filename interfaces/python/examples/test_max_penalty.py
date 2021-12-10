import lcqpow
import numpy as np


print("Preparing Maximum Penalty Test...")

# Setup data of first QP.
H = np.array([[2.0, 0.0], 
              [0.0, 2.0]])
g = np.array([-2.0, -2.0])
S1 = np.array([[1.0, 0.0]])
S2 = np.array([[0.0, 1.0]])

x0 = np.array([1.0, 1.0])
y0 = np.array([0.0, 0.0, 0.0, 0.0])

nV = 2 
nC = 0
nComp = 1
lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)

options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.QPOASES_DENSE)
options.setMaxRho(1.0)
lcqp.setOptions(options)

retVal = lcqp.loadLCQP(H=H, g=g, S1=S1, S2=S2, x0=x0, y0=y0)
if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal == lcqpow.ReturnValue.MAX_PENALTY_REACHED:
    print("Succeeded test (failed to solve LCQP due to maximal penalty).")