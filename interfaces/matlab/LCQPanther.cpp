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

#include <mex.h>

bool checkDimensionAndType(const mxArray* arr, int m, int n, const char* name) 
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
        sprintf(errorMsg, "Invalid dimension: %s (got %d x %d but expected %d x %d).\n", name, (int)mxGetM(arr), (int)mxGetN(arr), m, n);
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

void loadDenseArray(const mxArray* src, double* dest) 
{
    dest = (double*) mxGetPr( src );
}

void clearMemory(double* H, double* g, double* S1, double* S2, double* A, double* lbA, double* ubA, double* lb, double* ub)
{
    if (H != NULL) {
        delete[] H;
        H = NULL;
    }

    if (g != NULL) {
        delete[] g;
        g = NULL;
    } 

    if (S1 != NULL) {
        delete[] S1;
        S1 = NULL;
    }

    if (S2 != NULL) {
        delete[] S2;
        S2 = NULL;
    }

    if (lb != NULL) {
        delete[] lb;
        lb = NULL;
    }

    if (ub != NULL) {
        delete[] ub;
        ub = NULL;
    }

    if (A != NULL) {
        delete[] A;
        A = NULL;
    }

    if (lbA != NULL) {
        delete[] lbA;
        lbA = NULL;
    }

    if (ubA != NULL) {
        delete[] ubA;
        ubA = NULL;
    }
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    // Create a char for message handling
    char *errorMsg = (char*)malloc(100*sizeof(char));

    // Validate number of output arguments
    if (nlhs != 1) {
        sprintf(errorMsg, "Invalid number of output arguments (got %d but expected 1).\n", nlhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    // Validate number of input arguments
    if (nrhs < 4 || nrhs > 10) {
        sprintf(errorMsg, "Invalid number of input arguments (got %d but expected between 4 and 10).\n", nrhs);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    int nV = 0;
    int nComp = 0;
    int nC = 0;

    // Get number of optimization variables
    if (mxIsEmpty(prhs[0]) || !mxIsDouble(prhs[0])) {
        sprintf(errorMsg, "Invalid input argument: Hessian must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nV = (int) mxGetM(prhs[0]);
    }

    // Get number of complementarity constraints
    if (mxIsEmpty(prhs[2]) || !mxIsDouble(prhs[2])) {
        sprintf(errorMsg, "Invalid input argument: S1 must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nComp = (int) mxGetM(prhs[2]);
    }

    // Get number of linear constraints
    if (nrhs > 7 || ( nrhs == 7 && !mxIsStruct(prhs[6]))) {
        if (mxIsEmpty(prhs[4]) || !mxIsDouble(prhs[4])) {
            sprintf(errorMsg, "Invalid input argument: A must be a non-empty double matrix.\n");
            mexErrMsgTxt(errorMsg);
            free(errorMsg);
            return;
        } else {
            nC = (int) mxGetM(prhs[2]);
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
    LCQProblem lcqp(nV, nC, nComp);

    // Load data
    double* H = new double[nV*nV];
    double* g = new double[nV];
    double* S1 = new double[nComp*nV];
    double* S2 = new double[nComp*nV];
    double* lb = NULL;
    double* ub = NULL;
    double* A = NULL;
    double* lbA = NULL;
    double* ubA = NULL;
    // Options opts();
        
    loadDenseArray(prhs[0], H);
    loadDenseArray(prhs[1], g);
    loadDenseArray(prhs[2], S1);
    loadDenseArray(prhs[3], S2);

    if (nrhs == 6 || (nrhs == 7 && nC == 0)) {
        double* lb = new double[nV];
        double* ub = new double[nV];
        
        loadDenseArray(prhs[4], lb);
        loadDenseArray(prhs[5], ub);
    } else if (nrhs >= 7) {
        A = new double[nC*nV];
        lbA = new double[nC];
        ubA = new double[nC];

        loadDenseArray(prhs[4], A);
        loadDenseArray(prhs[5], lbA);
        loadDenseArray(prhs[6], ubA);

        if (nrhs >= 9) {
            double* lb = new double[nV];
            double* ub = new double[nV];
            
            loadDenseArray(prhs[7], lb);
            loadDenseArray(prhs[8], ub);
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

    // Load data into LCQP object
    lcqp.loadLCQP(H, g, S1, S2, A, lbA, ubA, lb, ub);

    // Run solver
    lcqpOASES::returnValue ret = lcqp.runSolver();
    if (ret != lcqpOASES::SUCCESSFUL_RETURN) {
        mexPrintf("Failed to solve LCQP (%d).\n", ret);
    } else {
        double* xOpt = new double[nV];
        lcqp.getPrimalSolution(xOpt);

        mexPrintf("Succeeded to solve LCQP. Obtaining solution vector.\n");
 
        plhs[0] = mxCreateDoubleMatrix(nV, 1, mxREAL);
        xOpt = mxGetPr(plhs[0]);
    }

    mexPrintf("Clearing memory.\n");

    // Clear all allocated memory
    clearMemory( H, g, S1, S2, A, lbA, ubA, lb, ub);    

    mexPrintf("Leaving mex function.\n");

    return;
}