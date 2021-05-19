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
#include <sys/stat.h>
#include <chrono>

bool PathExists(const std::string &s)
{
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
}

void WriteToFile(const double* const vec, int n, std::string path)
{
    // TODO; Delete file before open
    // std::ofstream myfile;
    // myfile.open(path);
    // for (int i = 0; i < n; i++)
    //     myfile << vec[i] << std::endl;
}

int main(int argc, char **argv) {

    // if (argc != 2) {
    //     printf("Wrong amount of arguments passed (not 2).");
    //     return -1;
    // }

    std::string inputdir = "examples/example_data/one_ivocp_example";

    if (!PathExists(inputdir)) {
        printf("Input directory does not exist.");
        return -1;
    }

    // Required files
    std::string H_file = inputdir + "/" + "H.txt";
    std::string g_file = inputdir + "/" + "g.txt";
    std::string lb_file = inputdir + "/" + "lb.txt";
    std::string ub_file = inputdir + "/" + "ub.txt";
    std::string S1_file = inputdir + "/" + "S1.txt";
    std::string S2_file = inputdir + "/" + "S2.txt";

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

    // Optional files
    std::string A_file = inputdir + "/" + "A.txt";
    std::string lbA_file = inputdir + "/" + "lbA.txt";
    std::string ubA_file = inputdir + "/" + "ubA.txt";
    std::string x0_file = inputdir + "/" + "x0.txt";
    std::string y0_file = inputdir + "/" + "y0.txt";

    const char* Af = 0;
    const char* lbAf = 0;
    const char* ubAf = 0;
    const char* x0f = 0;
    const char* y0f = 0;

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

    if (PathExists(x0_file)) {
        x0f = &x0_file[0];
    }

    if (PathExists(y0_file)) {
        y0f = &y0f[0];
    }

    lcqpOASES::LCQProblem lcqp( nV, nC, nComp );

    lcqpOASES::Options opts;
    opts.setPrintLevel( lcqpOASES::PrintLevel::INNER_LOOP_ITERATES );
    opts.setQPSolver( lcqpOASES::QPSolver::QPOASES_SPARSE);
    lcqp.setOptions( opts );

    // Run solver
	lcqpOASES::ReturnValue ret = lcqp.loadLCQP( &H_file[0], &g_file[0], &S1_file[0], &S2_file[0], Af, lbAf, ubAf, &lb_file[0], &ub_file[0], x0f, y0f );

    if (ret != lcqpOASES::SUCCESSFUL_RETURN)
    {
        printf("Failed to load LCQP.\n");
        return 1;
    }

    ret = lcqp.runSolver();

    if (ret != lcqpOASES::SUCCESSFUL_RETURN)
    {
        printf("Failed to solve LCQP.\n");
        return 1;
    }

    return 0;
}
