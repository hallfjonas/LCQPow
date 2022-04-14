/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
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


#include "LCQProblem.hpp"
#include <iostream>

using namespace LCQPow;

int main() {
    std::cout << "Preparing warm up problem with a simple inequality constraint...\n";

    /* Setup data of first QP. */
    double H[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double S1[1*2] = {1.0, 0.0};
    double S2[1*2] = {0.0, 1.0};
    double A[1*2] = {1.0, -1.0};
    double lbA[1] = { -0.5 };
    double ubA[1] = {  INFINITY };

    int nV = 2;
    int nC = 1;
    int nComp = 1;

    LCQProblem lcqp( nV, nC, nComp );

	Options options;
    options.setPrintLevel(PrintLevel::INNER_LOOP_ITERATES);
	lcqp.setOptions( options );

    // Solve first LCQP
	ReturnValue retVal = lcqp.loadLCQP( H, g, S1, S2, 0, 0, 0, 0, A, lbA, ubA );

    if (retVal != SUCCESSFUL_RETURN)
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

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to solve LCQP.\n");
        return 1;
    }

    // Get solutions
    double* xOpt = new double[2];
	double* yOpt = new double[nV + nC + 2*nComp];
	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );
	printf( "\nxOpt = [ %g, %g ];  yOpt = [ %g, %g, %g, %g, %g ]; \n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2],yOpt[3],yOpt[4] );

    delete[] xOpt; delete[] yOpt;

    return 0;
}
