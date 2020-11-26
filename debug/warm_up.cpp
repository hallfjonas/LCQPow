// Copyright 2020 Jonas Hall

#include <iostream>
#include <qpOASES.hpp>
USING_NAMESPACE_QPOASES

int main() {
    std::cout << "Preparing warm up problem...\n";

    /* Setup data of first QP. */
    real_t H[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    real_t g[2] = { -2.0, -2.0 };
    real_t lb[2] = { 0.0, 0.0 };
    real_t ub[2] = { 10000.0, 10000.0 };
    real_t C[2*2] = { 0.0, 1.0, 1.0, 0.0 };

    // Initial guess
    real_t x0[2] = { 1.0 + EPS, 1.0 };

    int_t nV = 2;

    LCQProblemB lcqp( nV );

	Options options;
    options.initialComplementarityPenalty = 1.5;
    options.complementarityPenaltyUpdate = 2;
	lcqp.setOptions( options );

	int_t nWSR = 10000000;

    // Solve first LCQP
	lcqp.init( H, g, lb, ub, C, nWSR, 0, x0 );

    // Get solutions
    real_t xOpt[2];
	real_t yOpt[2];
	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );
	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e ];  objVal = %e\n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],lcqp.getObjVal() );

    return 0;
}
