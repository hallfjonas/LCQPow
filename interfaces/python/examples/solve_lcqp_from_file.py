import lcqpow
import math


print("Preparing optimization problem loaded from file...")

# Setup data of first QP.
H_file = "example_data/one_ivocp_example/H.txt"
g_file = "example_data/one_ivocp_example/g.txt"
lb_file = "example_data/one_ivocp_example/lb.txt"
ub_file = "example_data/one_ivocp_example/ub.txt"
S1_file = "example_data/one_ivocp_example/S1.txt"
S2_file = "example_data/one_ivocp_example/S2.txt"

nV = 0 
with open(lb_file) as f:
   nV = sum([1 for line in f])

nComp = 0
with open(S1_file) as f:
   nComp = sum([1 for line in f])
nComp = math.floor(nComp/nV) 

A_file = "example_data/one_ivocp_example/A.txt"
lbA_file = "example_data/one_ivocp_example/lbA.txt"
ubA_file = "example_data/one_ivocp_example/ubA.txt"
x0_file = "example_data/one_ivocp_example/x0.txt"

nC = 0
with open(A_file) as f:
   nC = sum([1 for line in f])
nC = math.floor(nC/nV) 

lcqp = lcqpow.LCQProblem(nV=nV, nC=nC, nComp=nComp)
options = lcqpow.Options()
options.setPrintLevel(lcqpow.PrintLevel.INNER_LOOP_ITERATES)
options.setQPSolver(lcqpow.QPSolver.QPOASES_SPARSE)
lcqp.setOptions(options)

retVal = lcqp.loadLCQP(H_file=H_file, g_file=g_file, 
                       S1_file=S1_file, S2_file=S2_file, 
                       A_file=A_file, lbA_file=lbA_file, ubA_file=ubA_file, 
                       x0_file=x0_file)
if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to load LCQP.")

retVal = lcqp.runSolver()

if retVal != lcqpow.ReturnValue.SUCCESSFUL_RETURN:
    print("Failed to solve LCQP.")