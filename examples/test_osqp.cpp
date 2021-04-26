/*
 *	This file is part of lcqpOASES.
 *
 *	lcqpOASES -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	lcqpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	lcqpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with lcqpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include <iostream>
#include "LCQProblem.hpp"

using namespace lcqpOASES;

int main() {
    std::cout << "Preparing warm up problem...\n";

    /* Setup data of first QP. */
    int H_nnx = 3;
    double H_data[3] = { 2.0, 2.0 };
    int H_i[3] = {0, 1};
    int H_p[3] = {0, 1, 2};

    double g[2] = { -2.0, -2.0 };

    int A_nnx = 2;
    double A_data[2] = { 1.0, 1.0 };
    int A_i[2] = {0, 1};
    int A_p[3] = {0, 1, 2};

    double l[2] = {0.0, 0.0};
    double u[2] = {10000.0, 100000.0};

    int n = 2;
    int m = 2;

    // Exitflag
    int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    osqp_set_default_settings(settings);

    // Populate data
    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, H_nnx, H_data, H_i, H_p);
        data->q = g;
        data->A = csc_matrix(data->m, data->n, A_nnx, A_data, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    osqp_setup(&work, data, settings);

    exitflag = osqp_solve(work);

    printf("exitflag = %d\n", exitflag);

    OSQPSolution* sol(work->solution);
	printf( "\nxOpt = [ %g, %g ];  yOpt = [ %g, %g ]; \n\n",
			sol->x[0], sol->x[1], sol->y[0], sol->y[1]);

    return 0;
}
