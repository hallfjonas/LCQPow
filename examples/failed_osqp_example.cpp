/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2022 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include <iostream>
#include <fstream>
#include <LCQProblem.hpp>
#include <unistd.h>
#include <sys/stat.h>
#include <chrono>

bool PathExists(const std::string &s)
{
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
}

int main() {

    std::cout << "Preparing OCP loaded from file...\n";

    std::string inputdir = "examples/FailedOSQPExample";

    if (!PathExists(inputdir)) {
        printf("Input directory does not exist.");
        return -1;
    }

    // Required files
    std::string Q_file = inputdir + "/" + "Q.txt";
    std::string g_file = inputdir + "/" + "g.txt";
    std::string L_file = inputdir + "/" + "L.txt";
    std::string R_file = inputdir + "/" + "R.txt";
    std::string lb_L_file = inputdir + "/" + "lbL.txt";
    std::string lb_R_file = inputdir + "/" + "lbR.txt";
    std::string ub_L_file = inputdir + "/" + "ubL.txt";
    std::string ub_R_file = inputdir + "/" + "ubR.txt";

    int nV = 0;
    int nC = 0;
    int nComp = 0;

    // Count dimensions
    std::string line;

    std::ifstream Qfile(Q_file);
    while (std::getline(Qfile, line))
        nV++;
    nV = (int)sqrt(nV);

    std::ifstream Lfile(L_file);
    while (std::getline(Lfile, line))
        nComp++;

    nComp = nComp/nV;

    // Optional files
    std::string A_file = inputdir + "/" + "A.txt";
    std::string lbA_file = inputdir + "/" + "lbA.txt";
    std::string ubA_file = inputdir + "/" + "ubA.txt";

    const char* Af = 0;
    const char* lbAf = 0;
    const char* ubAf = 0;

    if (PathExists(A_file)) {
        Af = &A_file[0];

        std::ifstream Afile(A_file);
        while (std::getline(Afile, line))
            nC++;

        nC = nC/nV;
    }

    if (PathExists(lbA_file)) {
        lbAf = &lbA_file[0];
    }

    if (PathExists(ubA_file)) {
        ubAf = &ubA_file[0];
    }

    if (PathExists(lbA_file)) {
        lbAf = &lbA_file[0];
    }

    if (PathExists(ubA_file)) {
        ubAf = &ubA_file[0];
    }

    LCQPow::LCQProblem lcqp( nV, nC, nComp );

    LCQPow::Options options;
    options.setPrintLevel( LCQPow::PrintLevel::INNER_LOOP_ITERATES );
    
    // Switch to qpOASES to compare to succeeded convergence
    // options.setQPSolver( LCQPow::QPSolver::QPOASES_DENSE );
    options.setQPSolver( LCQPow::QPSolver::OSQP_SPARSE );
        
    // Turn off dynamic penalty updates
    options.setNDynamicPenalty(0);
    options.setMaxIterations(50);

    lcqp.setOptions( options );

    // Load data
	LCQPow::ReturnValue retVal = lcqp.loadLCQP( &Q_file[0], &g_file[0], &L_file[0], &R_file[0], &lb_L_file[0], NULL, &lb_R_file[0], NULL, Af, lbAf, ubAf);

    if (retVal != LCQPow::SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    // Must switch to sparse mode (if using a sparse solver)
    if (options.getQPSolver() >= LCQPow::QPOASES_SPARSE) {
        retVal = lcqp.switchToSparseMode( );

        if (retVal != LCQPow::SUCCESSFUL_RETURN)
        {
            printf("Failed to switch to sparse mode LCQP.\n");
            return 1;
        }
    }

    // Run the solver
    retVal = lcqp.runSolver();

    if (retVal != LCQPow::SUCCESSFUL_RETURN)
    {
        printf("Failed to solve LCQP (%d).\n", retVal);

        // Get QP solver exit flag
        if (retVal == LCQPow::ReturnValue::SUBPROBLEM_SOLVER_ERROR) {
            LCQPow::OutputStatistics stats;
            lcqp.getOutputStatistics(stats);

            int qpExitFlag = stats.getQPSolverExitFlag();
            printf("QP Solver exit flag: %d\n", qpExitFlag);
        }

        return 1;
    }

    return 0;
}
