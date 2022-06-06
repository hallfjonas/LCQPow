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
    std::cout << "Preparing dense warm up problem...\n";

    /* Setup data of first QP. */
    double Q[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double L[1*2] = {1.0, 0.0};
    double R[1*2] = {0.0, 1.0};

    double x0[2] = {1.0, 1.0};
    double y0[4] = {0.0, 0.0, 0.0, 0.0};

    int nV = 2;
    int nC = 0;
    int nComp = 1;

    // Load data
    LCQProblem lcqp( nV, nC, nComp );
	ReturnValue retVal = lcqp.loadLCQP( Q, g, L, R, 0, 0, 0, 0, 0, 0, 0, 0, 0, x0, y0 );
    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    // LCQPow options
	Options options;

    // Settings (OSQP)
    OSQPSettings* settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    osqp_set_default_settings(settings);
    settings->verbose = true;
    settings->polish = true;
    options.setOSQPOptions(settings);
    options.setQPSolver(QPSolver::OSQP_SPARSE);
    lcqp.setOptions( options );

    // Switch to sparse mode
    retVal = lcqp.switchToSparseMode( );
    if (retVal != LCQPow::SUCCESSFUL_RETURN)
    {
        printf("Failed to switch to sparse mode LCQP.\n");
        return 1;
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
