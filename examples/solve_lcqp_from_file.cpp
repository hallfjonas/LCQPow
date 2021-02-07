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
    std::string inputdir = "examples/example_data/";

    // Required files
    std::string H_file = inputdir + "H.txt";
    std::string g_file = inputdir + "g.txt";
    std::string lb_file = inputdir + "lb.txt";
    std::string ub_file = inputdir + "ub.txt";
    std::string S1_file = inputdir + "S1.txt";
    std::string S2_file = inputdir + "S2.txt";

    int nV = 0;
    int nC = 0;
    int nComp = 0;

    // Count dimensions
    std::string line;

    std::ifstream lbfile(lb_file);
    while (std::getline(lbfile, line))
        nV++;

    std::ifstream S1file(S1_file);
    while (std::getline(S1file, line))
        nComp++;

    nComp = nComp/nV;

    lcqpOASES::LCQProblem lcqp( nV, nC, nComp );
	lcqp.solve( &H_file[0], &g_file[0], &lb_file[0], &ub_file[0], &S1_file[0], &S2_file[0] );

    double* xOpt = new double[nV];
    double* yOpt = new double[nV];
    lcqp.getPrimalSolution( xOpt );
    lcqp.getDualSolution( yOpt );
}
