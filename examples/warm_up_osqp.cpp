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
    std::cout << "Preparing OSQP warm up problem...\n";

    /* Setup data of first QP. */
    double Q_data[2] = { 2.0, 2.0 };
    int Q_i[2] = {0, 1};
    int Q_p[3] = {0, 1, 2};

    double g[2] = { -2.0, -2.0 };

    double L_data[1] = { 1.0 };
    int L_i[1] = {0};
    int L_p[3] = {0, 1, 1};

    double R_data[1] = { 1.0 };
    int R_i[1] = {0};
    int R_p[3] = {0, 0, 1};

    int nV = 2;
    int nC = 0;
    int nComp = 1;

    csc* Q = LCQPow::Utilities::createCSC(nV, nV, Q_p[nV], Q_data, Q_i, Q_p);
    csc* L = LCQPow::Utilities::createCSC(nComp, nV, L_p[nV], L_data, L_i, L_p);
    csc* R = LCQPow::Utilities::createCSC(nComp, nV, R_p[nV], R_data, R_i, R_p);

    LCQProblem lcqp( nV, nC, nComp );

	Options options;
    options.setPrintLevel(PrintLevel::INNER_LOOP_ITERATES);
    options.setQPSolver(QPSolver::OSQP_SPARSE);
	lcqp.setOptions( options );

    // Solve first LCQP
	ReturnValue retVal = lcqp.loadLCQP(Q, g, L, R);

    free(Q); Q = NULL;
    free(L); L = NULL;
    free(R); R = NULL;

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    // Must switch to sparse mode (if using a sparse solver)
    if (options.getQPSolver() < LCQPow::QPOASES_SPARSE) {
        retVal = lcqp.switchToDenseMode( );

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

    int nDuals = lcqp.getNumberOfDuals();

    // Get solutions
    double* xOpt = new double[2];
	double* yOpt = new double[nDuals];
	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );

    if (nDuals == 2) {
        printf( "\nxOpt = [ %g, %g ];  yOpt = [ %g, %g ]; \n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1]);
    } else if (nDuals == 4) {
        printf( "\nxOpt = [ %g, %g ];  yOpt = [ %g, %g, %g, %g ]; \n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2],yOpt[3]);
    }

    delete[] xOpt; delete[] yOpt;

    return 0;
}
