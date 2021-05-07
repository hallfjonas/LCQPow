#include <mex.h>
#include "LCQProblem.hpp"

using lcqpOASES::LCQProblem;

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    int nV = 2;
    int nC = 0;
    int nComp = 1;

    mexPrintf("Creating LCQP objec.\n");

    LCQProblem lcqp( nV, nC, nComp );
    
    lcqp.~LCQProblem();

    mexPrintf("Leaving LCQParrot function.\n");
    return;
}