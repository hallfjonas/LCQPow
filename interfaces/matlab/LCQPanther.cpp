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

#include "LCQProblem.hpp"
using lcqpOASES::LCQProblem;
using lcqpOASES::Utilities;

#include <mex.h>

LCQProblem lcqp;

bool checkDimensionAndType(const mxArray* arr, size_t m, size_t n, const char* name)
{
    if (!mxIsDouble(arr)) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid type: %s must be of type double.\n", name);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return false;
    }

    if (mxGetM(arr) != m || mxGetN(arr) != n) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid dimension: %s (got %d x %d but expected %d x %d).\n", name, (int)mxGetM(arr), (int)mxGetN(arr), (int)m, (int)n);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return false;
    }

    return true;
}


bool checkStructType(const mxArray* arr, const char* name)
{
    if (!mxIsStruct(arr)) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid type: %s must be of type struct.\n", name);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return false;
    }

    return true;
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    // Validate number of output arguments
    if (nlhs != 0) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid number of output arguments (got %d but expected 1).\n", nlhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    // Validate number of input arguments
    if (nrhs < 4 || nrhs > 10) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid number of input arguments (got %d but expected between 4 and 10).\n", nrhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    size_t nV = 0;
    size_t nComp = 0;
    size_t nC = 0;

    // Get number of optimization variables
    if (mxIsEmpty(prhs[0]) || !mxIsDouble(prhs[0])) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid input argument: Hessian must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nV = mxGetM(prhs[0]);
    }

    // Get number of complementarity constraints
    if (mxIsEmpty(prhs[2]) || !mxIsDouble(prhs[2])) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid input argument: S1 must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nComp = mxGetM(prhs[2]);
    }

    // Get number of linear constraints
    if (nrhs > 7 || ( nrhs == 7 && !mxIsStruct(prhs[6]))) {
        if (mxIsEmpty(prhs[4]) || !mxIsDouble(prhs[4])) {
            char *errorMsg = (char*)malloc(100*sizeof(char));
            sprintf(errorMsg, "Invalid input argument: A must be a non-empty double matrix.\n");
            mexErrMsgTxt(errorMsg);
            free(errorMsg);
            return;
        } else {
            nC = mxGetM(prhs[2]);
        }
    }

    mexPrintf("Got LCQP of dimension (nV x nComp x nC) = (%d x %d x %d).\n", nV, nComp, nC);

    // Check all dimensions (except for params)
    if (!checkDimensionAndType(prhs[0], nV, nV, "H")) return;
    if (!checkDimensionAndType(prhs[1], nV, 1, "g")) return;
    if (!checkDimensionAndType(prhs[2], nComp, nV, "S1")) return;
    if (!checkDimensionAndType(prhs[3], nComp, nV, "S2")) return;

    if (nrhs == 6 && !checkDimensionAndType(prhs[4], nV, 1, "lb")) return;
    if (nrhs == 6 && !checkDimensionAndType(prhs[5], nV, 1, "ub")) return;

    if (nrhs == 7 && mxIsStruct(prhs[6]) && !checkDimensionAndType(prhs[4], nV, 1, "lb")) return;
    if (nrhs == 7 && mxIsStruct(prhs[6]) && !checkDimensionAndType(prhs[5], nV, 1, "ub")) return;

    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndType(prhs[4], nC, nV, "A")) return;
    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndType(prhs[5], nC, 1, "lbA")) return;
    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndType(prhs[6], nC, 1, "ubA")) return;

    if ((nrhs == 8 || nrhs == 9) && !checkDimensionAndType(prhs[4], nC, nV, "A")) return;
    if ((nrhs == 8 || nrhs == 9) && !checkDimensionAndType(prhs[5], nC, 1, "lbA")) return;
    if ((nrhs == 8 || nrhs == 9) && !checkDimensionAndType(prhs[6], nC, 1, "ubA")) return;

    if ((nrhs == 9 || nrhs == 10) && !checkDimensionAndType(prhs[7], nV, 1, "lb")) return;
    if ((nrhs == 9 || nrhs == 10) && !checkDimensionAndType(prhs[8], nV, 1, "ub")) return;

    // Check structs
    if (nrhs == 5 && !checkStructType(prhs[4], "params")) return;
    if (nrhs == 7 && !mxIsDouble(prhs[6]) && !checkStructType(prhs[6], "params")) return;
    if (nrhs == 8 && !checkStructType(prhs[7], "params")) return;
    if (nrhs == 10 && !checkStructType(prhs[9], "params")) return;

    // Sparsity checks
    if (mxIsSparse(prhs[0]) || (nC > 0 && mxIsSparse(prhs[4]))) {
        mexErrMsgTxt("Sparsity is not yet suported.\n");
        return;
    }

    // Create LCQP object
    LCQProblem lcqp((int)nV, (int)nC, (int)nComp);

    // Load data
    double* H = NULL;
    double* g = NULL;
    double* S1 = NULL;
    double* S2 = NULL;
    double* lb = NULL;
    double* ub = NULL;
    double* A = NULL;
    double* lbA = NULL;
    double* ubA = NULL;
    // Options opts();

    H = (double*) mxGetPr( prhs[0] );
    g = (double*) mxGetPr( prhs[1] );
    S1 = (double*) mxGetPr( prhs[2] );
    S2 = (double*) mxGetPr( prhs[3] );

    if (nrhs == 6 || (nrhs == 7 && nC == 0)) {
        lb = (double*) mxGetPr( prhs[4] );
        ub = (double*) mxGetPr( prhs[5] );
    } else if (nrhs >= 7) {
        A = (double*) mxGetPr( prhs[4] );
        lbA = (double*) mxGetPr( prhs[5] );
        ubA = (double*) mxGetPr( prhs[6] );

        if (nrhs >= 9) {
            lb = (double*) mxGetPr( prhs[7] );
            ub = (double*) mxGetPr( prhs[8] );
        }
    }

    // Load settings
    int structIdx = -1;
    if (nrhs == 5 && checkStructType(prhs[4], "params")) { structIdx = 4; }
    if (nrhs == 7 && !mxIsDouble(prhs[6]) && checkStructType(prhs[6], "params"))  { structIdx = 6; }
    if (nrhs == 8 && checkStructType(prhs[7], "params")) { structIdx = 7; }
    if (nrhs == 10 && checkStructType(prhs[9], "params")) { structIdx = 9; }

    if (structIdx != -1) {
        mexErrMsgTxt("Passing the parameter field is not yet supported.\n");
        return;
    }

    Utilities::printMatrix(H, nV, nV, "H");
    Utilities::printMatrix(g, 1, nV, "H");
    Utilities::printMatrix(S1, nComp, nV, "H");
    Utilities::printMatrix(S2, nComp, nV, "H");

    // double* xOpt;

    // Load data into LCQP object
    lcqp.loadLCQP(H, g, S1, S2, A, lbA, ubA, lb, ub);

    // Run solver
    lcqpOASES::returnValue ret = lcqp.runSolver();
    if (ret != lcqpOASES::SUCCESSFUL_RETURN) {
        mexPrintf("Failed to solve LCQP (%d).\n", ret);
    } else {
        mexPrintf("Succeeded to solve LCQP. Obtaining solution vector.\n");

        // double* xOpt_tmp = new double[nV];
        // lcqp.getPrimalSolution(xOpt_tmp);


        // Utilities::printMatrix(xOpt_tmp, 1, nV, "xOpt");

        // delete[] xOpt_tmp;

        // plhs[0] = mxCreateDoubleMatrix(nV, 1, mxREAL);
        // if (plhs[0] == NULL) {
        //     mexPrintf("Failed to allocate output.\n");
        //     return;
        // }

        // xOpt = (double*) mxGetPr(plhs[0]);

        // for (int i = 0; i < nV; i++)
        //     xOpt[i] = xOpt_tmp[i];
    }

    // Destroy LCQProblem object
    lcqp.~LCQProblem();

    return;
}