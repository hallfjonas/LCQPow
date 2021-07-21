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


#include <iostream>
#include "LCQProblem.hpp"

using namespace LCQPow;

int main() {
    std::cout << "Preparing dense warm up problem...\n";

    /* Setup data of first QP. */
    double H[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };

    // 0 <= x _|_ y >= 0
    // 0 <= x _|_ 0.5 - x >= 0
    double S1[2*2] = {1.0, 0.0, 1.0, 0.0};
    double S2[2*2] = {0.0, 1.0, -1.0, 0.0};
    double lbS1[2] = {0.0, 0.0};
    double lbS2[2] = {0.0, -0.5};

    double x0[2] = {0.0, 0.0};

    int nV = 2;
    int nC = 0;
    int nComp = 2;

    LCQProblem lcqp( nV, nC, nComp );

	Options options;
    options.setPrintLevel(PrintLevel::INNER_LOOP_ITERATES);
    options.setQPSolver(QPSolver::QPOASES_DENSE);
	lcqp.setOptions( options );

    // Solve first LCQP
	ReturnValue retVal = lcqp.loadLCQP( H, g, S1, S2, lbS1, 0, lbS2, 0, 0, 0, 0, 0, 0, x0);

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    retVal = lcqp.runSolver();

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to solve LCQP.\n");
        return 1;
    }

    // Solve another LCQP
    x0[0] = 0.0;
    x0[1] = 3000;
    options.setSolveZeroPenaltyFirst(false);
    options.setInitialPenaltyParameter(10.0);
	lcqp.setOptions( options );
	retVal = lcqp.loadLCQP( H, g, S1, S2, lbS1, 0, lbS2, 0, 0, 0, 0, 0, 0, x0);

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    retVal = lcqp.runSolver();

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to solve LCQP.\n");
        return 1;
    }

    // Get solutions
    double* xOpt = new double[2];
	double* yOpt = new double[nV + nC + 2*nComp];
    LCQPow::OutputStatistics stats;
	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );
    lcqp.getOutputStatistics( stats );
	printf( "\nxOpt = [ %g, %g ];  yOpt = [ %g, %g, %g, %g ]; i = %d; k = %d; rho = %g; WSR = %d \n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2],yOpt[3],
            stats.getIterTotal(), stats.getIterOuter(), stats.getRhoOpt(), stats.getSubproblemIter() );

    // Clean Up
    delete[] xOpt; delete[] yOpt;

    return 0;
}
