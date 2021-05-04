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

#include <mex.h>


bool checkDimensionAndType(const mxArray* arr, int m, int n, const char* name) 
{
    if (!mxIsDouble(arr)) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid type: %s must be of type double", name);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return false;
    }

    if (mxGetM(arr) != m || mxGetN(arr) != n) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid dimension: %s (got %d x %d but expected %d x %d)", name, mxGetM(arr), mxGetN(arr), m, n);
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
        sprintf(errorMsg, "Invalid type: %s must be of type struct", name);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return false;
    }

    return true;
}


void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    // Create a char for message handling
    char *errorMsg = (char*)malloc(100*sizeof(char));

    // Validate number of output arguments
    if (nlhs != 1) {
        sprintf(errorMsg, "Invalid number of output arguments (got %d but expected 1)", nlhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    // Validate number of input arguments
    if (nrhs < 4 || nrhs > 10) {
        sprintf(errorMsg, "Invalid number of input arguments (got %d but expected between 4 and 10).", nrhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    int nV = 0;
    int nComp = 0;
    int nC = 0;

    // Get number of optimization variables
    if (mxIsEmpty(prhs[0]) || !mxIsDouble(prhs[0])) {
        sprintf(errorMsg, "Invalid input argument: Hessian must be a non-empty double matrix.", nrhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nV = (int) mxGetM(prhs[0]);
    }

    // Get number of complementarity constraints
    if (mxIsEmpty(prhs[2]) || !mxIsDouble(prhs[2])) {
        sprintf(errorMsg, "Invalid input argument: S1 must be a non-empty double matrix.", nrhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nComp = (int) mxGetM(prhs[2]);
    }

    // Get number of linear constraints
    if (nrhs > 7 || ( nrhs == 7 && !mxIsStruct(prhs[6]))) {
        if (mxIsEmpty(prhs[4]) || !mxIsDouble(prhs[4])) {
            sprintf(errorMsg, "Invalid input argument: A must be a non-empty double matrix.", nrhs);
            mexErrMsgTxt(errorMsg);
            free(errorMsg);
            return;
    } else {
        nC = (int) mxGetM(prhs[2]);
    }

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

    // To know we ran successfully.
    mexPrintf("Leaving mex function.");
    return;
}