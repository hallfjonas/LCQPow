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
#include <fstream>
#include <LCQProblem.hpp>
#include <unistd.h>

int main() {
    char* inputdir = "examples/LCQP/example_data/";

    // Required files
    const char* H_file = strcat(inputdir, "H.txt");
    const char* g_file = strcat(inputdir, "g.txt");
    const char* lb_file = strcat(inputdir, "lb.txt");
    const char* ub_file = strcat(inputdir, "ub.txt");
    const char* S1_file = strcat(inputdir, "S1.txt");
    const char* S2_file = strcat(inputdir, "S2.txt");

    // Contraints (optional files, but if A exists then all are required)
    const char* A_file = strcat(inputdir, "A.txt");
    const char* lbA_file = strcat(inputdir, "lbA.txt");
    const char* ubA_file = strcat(inputdir, "ubA.txt");

    // Initial primal and dual guess (optional)
    const char* x0_file = strcat(inputdir, "x0.txt");
    const char* y0_file = strcat(inputdir, "y0.txt");

    int nV = 0;
    int nC = 0;
    int nComp = 0;

    // Count dimensions
    std::string line;

    std::ifstream lbfile(lb_file);
    while (std::getline(lbfile, line))
        nV++;

    std::ifstream lbAfile(lbA_file);
    while (std::getline(lbAfile, line))
        nC++;

    std::ifstream S1file(S1_file);
    while (std::getline(S1file, line))
        nComp++;

    lcqpOASES::LCQProblem lcqp( nV, nC, nComp );
	lcqp.solve( H_file, g_file, lb_file, ub_file, S1_file, S2_file );

    double* xOpt = new double[nV];
    double* yOpt = new double[nV];
    lcqp.getPrimalSolution( xOpt );
    lcqp.getDualSolution( yOpt );
}
