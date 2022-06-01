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
#include "LCQProblem.hpp"

using namespace LCQPow;

int main() {
    std::cout << "Preparing unit circle optimization problem...\n";

    // Set dimensions
    int N = 100;
    int nV = 2 + 2*N;
    int nC = N+1;
    int nComp = N;

    // Reference
    double x_ref[2] = {0.5, -0.6};

    // Set up LCQP object
    LCQProblem lcqp( nV, nC, nComp );
	Options options;
    options.setPrintLevel(PrintLevel::SUBPROBLEM_SOLVER_ITERATES);
    options.setQPSolver(QPSolver::QPOASES_SPARSE);
    // options.setStationarityTolerance( 10e-5 );
    options.setStoreSteps(true);
    lcqp.setOptions( options );

    // Allocate vectors
    double* Q = new double[nV*nV]();
    double* g = new double[nV]();
    double* L = new double[nComp*nV]();
    double* R = new double[nComp*nV]();
    double* A = new double[nC*nV]();
    double* lbA = new double[nC]();
    double* ubA = new double[nC]();
    double* x0 = new double[nV]();
    x0[0] = x_ref[0];
    x0[1] = x_ref[1];

    // Assign problem data
    Q[0] = 17; Q[1*nV + 1] = 17;
    Q[0*nV + 1] = -15;
    Q[1*nV + 0] = -15;

    // Regularization on Q
    for (int i = 2; i < nV; i++)
        Q[i*nV + i] = 5e-12;

    // Objective linear term
    double Qx[2*2] = {17, -15, -15, 17};
    LCQPow::Utilities::MatrixMultiplication(Qx, x_ref, g, 2, 2, 1);
    LCQPow::Utilities::WeightedVectorAdd(-1, g, 0, g, g, 2);

    // Constraints
    for (int i = 0; i < N; i++) {
        // Equality constraint ([cos sin]'*x + lambda = 1)
        A[i*nV + 0] = cos((2*M_PI*i)/N);
        A[i*nV + 1] = sin((2*M_PI*i)/N);
        A[i*nV + 2 + 2*i] = 1;

        // Convex combination constraint (sum theta = 1)
        A[N*nV + 3 + 2*i] = 1;

        // Complementarity constraints (lamba*theta = 0)
        L[i*nV + 2 + 2*i] = 1;
        R[i*nV + 3 + 2*i] = 1;

        // constraint bounds
        lbA[i] = 1;
        ubA[i] = 1;

        x0[2*i + 2] = 1;
        x0[2*i + 3] = 1;
    }

    // Constraint bound for last constraint
    lbA[N] = 1;
    ubA[N] = 1;

    // Load Data
	ReturnValue retVal = lcqp.loadLCQP( Q, g, L, R, 0, 0, 0, 0,  A, lbA, ubA, 0, 0, x0 );

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    // Must switch to sparse mode (if using a sparse solver)
    if (options.getQPSolver() >= LCQPow::QPOASES_SPARSE) {
        retVal = lcqp.switchToSparseMode( );

        if (retVal != SUCCESSFUL_RETURN)
        {
            printf("Failed to switch to sparse mode LCQP.\n");
            return 1;
        }
    }

    // Run solver
    retVal = lcqp.runSolver();

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to solve LCQP.\n");
        return 1;
    }

    // Get solutions
    double* xOpt = new double[nV];
	double* yOpt = new double[nV + nC + 2*nComp];
    LCQPow::OutputStatistics stats;
	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );
    lcqp.getOutputStatistics( stats );

    // Print coordinates of solution and output statistics
	printf( "\nxOpt = [ %g, %g ];  i = %d; k = %d; rho = %g; WSR = %d \n\n",
			xOpt[0],xOpt[1],
            stats.getIterTotal(), stats.getIterOuter(), stats.getRhoOpt(), stats.getSubproblemIter() );

    // Print a reference to the global and local solutions
    printf("For reference: Global solution is at:  [ %g, %g ]\n", 0.1811, -0.9835);
    printf("               Another local solution: [ %g, %g ]\n",  0.9764, -0.2183);

    int* innerItersTMP = stats.getInnerIters( );
    int* subproblemItersTMP = stats.getSubproblemIters( );
    double* stepLengthTMP = stats.getStepLength( );
    double* stepSizeTMP = stats.getStepSize( );
    double* statValsTMP = stats.getStatVals( );
    double* objValsTMP = stats.getObjVals( );
    double* phiValsTMP = stats.getPhiVals( );
    double* meritValsTMP = stats.getMeritVals( );

    printf(" %6s |", "inner");
    printf(" %6s |", "subpr.");
    printf(" %10s |", "alpha");
    printf(" %10s |", "norm pk");
    printf(" %10s |", "statVals");
    printf(" %10s |", "objVals");
    printf(" %10s |", "phiVals");
    printf(" %10s\n", "meritVals");

    for (int i = 0; i < stats.getIterTotal()+1; i++) {
        printf(" %6d |", innerItersTMP[i]);
        printf(" %6d |", subproblemItersTMP[i]);
        printf(" %10.4g |", stepLengthTMP[i]);
        printf(" %10.4g |", stepSizeTMP[i]);
        printf(" %10.4g |", statValsTMP[i]);
        printf(" %10.4g |", objValsTMP[i]);
        printf(" %10.4g |", phiValsTMP[i]);
        printf(" %10.4g\n", meritValsTMP[i]);
    }

    // Clean Up
    delete[] xOpt; delete[] yOpt;
    delete[] Q; delete[] g; delete[] L; delete[] R; delete[] A;
    delete[] lbA; delete[] ubA; delete[] x0;

    delete[] innerItersTMP;
    delete[] subproblemItersTMP;
    delete[] stepLengthTMP;
    delete[] stepSizeTMP;
    delete[] statValsTMP;
    delete[] objValsTMP;
    delete[] phiValsTMP;
    delete[] meritValsTMP;

    return 0;
}
