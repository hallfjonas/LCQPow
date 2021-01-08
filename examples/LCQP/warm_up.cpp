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
    real_t S1[1*2] = {1.0, 0.0};
    real_t S2[1*2] = {0.0, 1.0};

    real_t A[1*2] = { 1.0, 0.0 };
    real_t lbA[1] = { -10.0 };
    real_t ubA[1] = { 100.0};

    // Initial guess
    real_t x0[2] = { 1.0 + EPS, 1.0 };

    int_t nV = 2;
    int_t nC = 1;
    int_t nComp = 1;

    LCQProblem lcqp( nV, nC, nComp );

	Options options;
    options.initialComplementarityPenalty = 1.5;
    options.complementarityPenaltyUpdate = 2;
	lcqp.setOptions( options );

	int_t nWSR = 10000000;

    // Solve first LCQP
	returnValue retVal = lcqp.solve( H, g, A, lb, ub, lbA, ubA, S1, S2, nWSR, 0, x0 );

    if (retVal != SUCCESSFUL_RETURN)
    {
        printf("Termination ended without success.\n");
        return 1;
    }

    // Get solutions
    real_t xOpt[2];
	real_t yOpt[2];
	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );
	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e ];  objVal = %e\n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],lcqp.getObjVal() );

    return 0;
}
