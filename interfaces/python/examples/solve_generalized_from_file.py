import lcqpow
import math


print("Preparing optimization problem loaded from file...")

# Setup data of first QP.
Q_file = "example_data/generalized_constraints/Q.txt"
g_file = "example_data/generalized_constraints/g.txt"
lb_file = "example_data/generalized_constraints/lb.txt"
ub_file = "example_data/generalized_constraints/ub.txt"
L_file = "example_data/generalized_constraints/L.txt"
R_file = "example_data/generalized_constraints/R.txt"
lbL_file = "example_data/generalized_constraints/lbL.txt"
lbR_file = "example_data/generalized_constraints/lbR.txt"
ubL_file = "example_data/generalized_constraints/ubL.txt"
ubR_file = "example_data/generalized_constraints/ubR.txt"

nV = 0 
with open(lb_file) as f:
   nV = sum([1 for line in f])

nComp = 0
with open(L_file) as f:
   nComp = sum([1 for line in f])
nComp = math.floor(nComp/nV) 

A_file = "example_data/generalized_constraints/A.txt"
lbA_file = "example_data/generalized_constraints/lbA.txt"
ubA_file = "example_data/generalized_constraints/ubA.txt"
x0_file = "example_data/generalized_constraints/x0.txt"

nC = 0
with open(A_file) as f:
   nC = sum([1 for line in f])
nC = math.floor(nC/nV) 

lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)
options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.QPOASES_SPARSE)
lcqp.setOptions(options)

retVal = lcqp.loadLCQP(Q_file=Q_file, g_file=g_file, 
                       L_file=L_file, R_file=R_file, 
                       lbL_file=lbL_file, ubL_file=ubL_file, 
                       lbR_file=lbR_file, ubR_file=ubR_file, 
                       A_file=A_file, lbA_file=lbA_file, ubA_file=ubA_file, 
                       x0_file=x0_file)
if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to solve LCQP.")