#include <mex.h>
#include "LCQProblem.hpp"

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    int nV = 2;
    int nC = 0;
    int nComp = 1;

    mexPrintf("Creating LCQP objec.\n");

    lcqpOASES::LCQProblem lcqp( nV, nC, nComp );

    mexPrintf("Leaving LCQParrot function.\n");
    return;
}