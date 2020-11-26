// Copyright 2020 Jonas Hall

#include <iostream>
#include <fstream>
#include <qpOASES.hpp>
#include <unistd.h>

USING_NAMESPACE_QPOASES

int main() {
    std::string inputdir = "examples/LCQP/example_data/";

    // Required files
    std::string H_file = inputdir + "H.txt";
    std::string g_file = inputdir + "g.txt";
    std::string lb_file = inputdir + "lb.txt";
    std::string ub_file = inputdir + "ub.txt";
    std::string C_file = inputdir + "C.txt";

    // Contraints (optional files, but if A exists then all are required)
    std::string A_file = inputdir + "A.txt";
    std::string lbA_file = inputdir + "lbA.txt";
    std::string ubA_file = inputdir + "ubA.txt";

    // Initial primal and dual guess (optional)
    std::string x0_file = inputdir + "x0.txt";
    std::string y0_file = inputdir + "y0.txt";

    int_t nV = 0;
    int_t nC = 0;

    std::string line;

    std::ifstream lbfile(lb_file);
    while (std::getline(lbfile, line))
        nV++;

    std::ifstream lbAfile(lbA_file);
    while (std::getline(lbAfile, line))
        nC++;

    // Read x0 value
    real_t* x0 = new real_t[nV];
    if (access( x0_file.c_str(), F_OK ) != -1 ) {
        readFromFile( x0, nV, &x0_file[0] );
    } else {
        for (int i = 0; i < nV; i++)
            x0[i] = 0;
    }

    // Read y0 value
    real_t* y0 = new real_t[nV + nC];
    if (access( y0_file.c_str(), F_OK ) != -1 ) {
        readFromFile( y0, nV + nC, &y0_file[0] );
    } else {
        for (int i = 0; i < nV + nC; i++)
            y0[i] = 0;
    }

    Options options;
    options.setToDefault();
    options.initialComplementarityPenalty = 1.0;
    options.complementarityPenaltyUpdate = 2.0;

    int_t nWSR = 10000000;
    real_t* cputime = 0;

    if (nC == 0) {
        LCQProblemB lcqp( nV );
	    lcqp.setOptions( options );
        lcqp.init( &H_file[0], &g_file[0], &lb_file[0], &ub_file[0], &C_file[0], nWSR, cputime, x0, y0 );

        real_t* xOpt = new real_t[nV];
        real_t* yOpt = new real_t[nV];
        lcqp.getPrimalSolution( xOpt );
        lcqp.getDualSolution( yOpt );
    } else {
        LCQProblem lcqp( nV, nC );
	    lcqp.setOptions( options );
        lcqp.init( &H_file[0], &g_file[0], &A_file[0], &lb_file[0], &ub_file[0], &lbA_file[0], &ubA_file[0], &C_file[0], nWSR, cputime, x0, y0 );

        real_t* xOpt = new real_t[nV];
        real_t* yOpt = new real_t[nV];
        lcqp.getPrimalSolution( xOpt );
        lcqp.getDualSolution( yOpt );
    }
}
