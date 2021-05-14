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
using lcqpOASES::Options;

#include <mex.h>
#include <chrono>

LCQProblem lcqp;

bool checkDimensionAndTypeDouble(const mxArray* arr, size_t m, size_t n, const char* name)
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

bool checkDimensionAndTypeBool(const mxArray* arr, size_t m, size_t n, const char* name)
{
    if (!mxIsLogical(arr)) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid type: %s must be of type logical.\n", name);
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

bool checkTypeStruct(const mxArray* arr, const char* name)
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

void colMajorToRowMajor(double* col_maj, double* row_maj, int m, int n)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            row_maj[i*n + j] = col_maj[j*m + i];
}

void printOptions( Options options ) {
    mexPrintf(" \n Using LCQPanther Options: \n");
    mexPrintf("          rho0: %g \n", options.getInitialPenaltyParameter());
    mexPrintf("          beta: %g \n", options.getPenaltyUpdateFactor());
    mexPrintf("     compl tol: %g \n", options.getComplementarityTolerance());
    mexPrintf("     stati tol: %g \n", options.getStationarityTolerance());
    mexPrintf("      max iter: %d \n", options.getMaxIterations());
    mexPrintf("solve zero pen: %d \n", (int) options.getSolveZeroPenaltyFirst());
    mexPrintf("          prnt: %d \n\n", options.getPrintLevel());
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    // Validate number of output arguments
    int nlhs_min = 1; int nlhs_max = 3;
    if (nlhs < nlhs_min || nlhs > nlhs_max) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid number of output arguments (got %d but expected between %d and %d).\n", nlhs, nlhs_min, nlhs_max);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    }

    // Validate number of input arguments
    int nrhs_min = 1; int nrhs_max = 10;
    if (nrhs < nrhs_min || nrhs > nrhs_max) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid number of input arguments (got %d but expected between %d and %d).\n", nrhs, nrhs_min, nrhs_max);
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
            nC = mxGetM(prhs[4]);
        }
    }

    // Debug statement
    // mexPrintf("Got LCQP of dimension (nV x nComp x nC) = (%d x %d x %d).\n", nV, nComp, nC);

    // Check all dimensions (except for params)
    if (!checkDimensionAndTypeDouble(prhs[0], nV, nV, "H")) return;
    if (!checkDimensionAndTypeDouble(prhs[1], nV, 1, "g")) return;
    if (!checkDimensionAndTypeDouble(prhs[2], nComp, nV, "S1")) return;
    if (!checkDimensionAndTypeDouble(prhs[3], nComp, nV, "S2")) return;

    if (nrhs == 6 && !checkDimensionAndTypeDouble(prhs[4], nV, 1, "lb")) return;
    if (nrhs == 6 && !checkDimensionAndTypeDouble(prhs[5], nV, 1, "ub")) return;

    if (nrhs == 7 && mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[4], nV, 1, "lb")) return;
    if (nrhs == 7 && mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[5], nV, 1, "ub")) return;

    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[4], nC, nV, "A")) return;
    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[5], nC, 1, "lbA")) return;
    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[6], nC, 1, "ubA")) return;

    if ((nrhs == 8 || nrhs == 9) && !checkDimensionAndTypeDouble(prhs[4], nC, nV, "A")) return;
    if ((nrhs == 8 || nrhs == 9) && !checkDimensionAndTypeDouble(prhs[5], nC, 1, "lbA")) return;
    if ((nrhs == 8 || nrhs == 9) && !checkDimensionAndTypeDouble(prhs[6], nC, 1, "ubA")) return;

    if ((nrhs == 9 || nrhs == 10) && !checkDimensionAndTypeDouble(prhs[7], nV, 1, "lb")) return;
    if ((nrhs == 9 || nrhs == 10) && !checkDimensionAndTypeDouble(prhs[8], nV, 1, "ub")) return;

    // Check structs
    if (nrhs == 5 && !checkTypeStruct(prhs[4], "params")) return;
    if (nrhs == 7 && !mxIsDouble(prhs[6]) && !checkTypeStruct(prhs[6], "params")) return;
    if (nrhs == 8 && !checkTypeStruct(prhs[7], "params")) return;
    if (nrhs == 10 && !checkTypeStruct(prhs[9], "params")) return;

    // Sparsity checks
    if (mxIsSparse(prhs[0]) || (nC > 0 && mxIsSparse(prhs[4]))) {
        mexErrMsgTxt("Sparsity is not yet suported.\n");
        return;
    }

    // Create LCQP and options objects
    LCQProblem lcqp((int)nV, (int)nC, (int)nComp);
    Options options;

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
    double* x0 = NULL;
    double* y0 = NULL;

    // Temporary matrices to keep the column major formats
    double* S1_col = NULL;
    double* S2_col = NULL;
    double* A_col = NULL;

    H = (double*) mxGetPr( prhs[0] );
    g = (double*) mxGetPr( prhs[1] );
    S1_col = (double*) mxGetPr( prhs[2] );
    S2_col = (double*) mxGetPr( prhs[3] );

    if (nrhs == 6 || (nrhs == 7 && nC == 0)) {
        lb = (double*) mxGetPr( prhs[4] );
        ub = (double*) mxGetPr( prhs[5] );
    } else if (nrhs >= 7) {
        A_col = (double*) mxGetPr( prhs[4] );
        lbA = (double*) mxGetPr( prhs[5] );
        ubA = (double*) mxGetPr( prhs[6] );

        if (nrhs >= 9) {
            lb = (double*) mxGetPr( prhs[7] );
            ub = (double*) mxGetPr( prhs[8] );
        }
    }

    // MATLAB stores in column major format (switch to row major)
    if (S1_col != NULL && nComp > 1 && nV > 1) {
        S1 = new double[nComp*nV];
        colMajorToRowMajor(S1_col, S1, nComp, nV);
    }
    if (S2_col != NULL && nComp > 1 && nV > 1) {
        S2 = new double[nComp*nV];
        colMajorToRowMajor(S2_col, S2, nComp, nV);
    }
    if (A_col != NULL && nC > 1 && nV > 1) {
        A = new double[nC*nV];
        colMajorToRowMajor(A_col, A, nC, nV);
    }

    // Load settings
    int structIdx = -1;
    if (nrhs == 5 && checkTypeStruct(prhs[4], "params")) { structIdx = 4; }
    if (nrhs == 7 && !mxIsDouble(prhs[6]) && checkTypeStruct(prhs[6], "params"))  { structIdx = 6; }
    if (nrhs == 8 && checkTypeStruct(prhs[7], "params")) { structIdx = 7; }
    if (nrhs == 10 && checkTypeStruct(prhs[9], "params")) { structIdx = 9; }

    if (structIdx != -1) {

        const char* params_fieldnames[] = {
            "x0",
            "stationarityTolerance",
            "complementarityTolerance",
            "initialPenaltyParameter",
            "penaltyUpdateFactor",
            "solveZeroPenaltyFirst",
            "maxIterations",
            "printLevel"
        };

        for (auto name : params_fieldnames) {
            mxArray* field = mxGetField(prhs[structIdx], 0, name);
            double* fld_ptr;
            bool* fld_ptr_bool;

            if (field == NULL)
                continue;

            if ( name == "stationarityTolerance") {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.stationarityTolerance")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setStationarityTolerance( fld_ptr[0] );
                continue;
            }

            if ( name == "complementarityTolerance") {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.complementarityTolerance")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setComplementarityTolerance( fld_ptr[0] );
                continue;
            }

            if ( name == "initialPenaltyParameter") {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.initialPenaltyParameter")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setInitialPenaltyParameter( fld_ptr[0] );
                continue;
            }

            if ( name == "penaltyUpdateFactor") {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.penaltyUpdateFactor")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setPenaltyUpdateFactor( fld_ptr[0] );
                continue;
            }

            if ( name == "solveZeroPenaltyFirst") {
                if (!checkDimensionAndTypeBool(field, 1, 1, "params.solveZeroPenaltyFirst")) return;

                fld_ptr_bool = (bool*) mxGetPr(field);
                options.setSolveZeroPenaltyFirst( fld_ptr[0] );
                continue;
            }

            if ( name == "maxIterations") {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.maxIterations")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setMaxIterations( (int)fld_ptr[0] );
                continue;
            }

            if ( name == "printLevel") {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.printLevel")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setPrintLevel( (int)fld_ptr[0] );
                continue;
            }

            if ( name == "x0") {
                if (!checkDimensionAndTypeDouble(field, nV, 1, "params.x0")) return;

                x0 = (double*) mxGetPr(field);
                continue;
            }
        }

        lcqp.setOptions( options );
    }

    // printOptions( options );

    // Load data into LCQP object
    lcqp.loadLCQP(H, g, S1, S2, A, lbA, ubA, lb, ub, x0);

    // Clear A, S1, S2
    if (A != 0)
        delete[] A;

    if (S1 != 0)
        delete[] S1;

    if (S2 != 0)
        delete[] S2;


    // Start time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Run solver
    lcqpOASES::ReturnValue ret = lcqp.runSolver();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsed_secs = (end - begin).count()/1000.0/1000.0/1000.0;

    if (ret != lcqpOASES::SUCCESSFUL_RETURN) {
        mexPrintf("Failed to solve LCQP (error code: %d).\n", ret);
    } else {
        // Succeeded to solve LCQP. Obtaining solution vector
        // 1) Primal solution vector
        plhs[0] = mxCreateDoubleMatrix(nV, 1, mxREAL);
        if (plhs[0] == NULL) {
            mexPrintf("Failed to allocate output of primal solution vector.\n");
            return;
        }

        // Point to the output object
        double* xOpt = (double*) mxGetPr(plhs[0]);
        lcqp.getPrimalSolution(xOpt);

        // 2) Dual solution vector
        if (nlhs > 1) {
            int nDuals = lcqp.getNumerOfDuals();

            if (nDuals <= 0) {
                mexPrintf("Failed to receive number of dual variables.\n");
                return;
            }

            plhs[1] = mxCreateDoubleMatrix(nDuals, 1, mxREAL);
            if (plhs[1] == NULL) {
                mexPrintf("Failed to allocate output of dual solution vector.\n");
                return;
            }

            // Pointer to the output object
            double* yOpt = (double*) mxGetPr(plhs[1]);
            lcqp.getDualSolution(yOpt);
        }

        // 3) Output statistics
        if (nlhs > 2) {
            // assign fieldnames
            int numberStatOutputs = 5;

            const char* fieldnames[] = {"iters_total", "iters_outer", "iters_subproblem", "rho_opt", "elapsed_time"};

            // Allocate memory
            plhs[2] = mxCreateStructMatrix(1,1,numberStatOutputs, fieldnames);

            // assert that allocation went ok
            if (plhs[2] == NULL) {
                mexPrintf("Failed to allocate output of statistics struct.\n");
                return;
            }

            // Get the statistics
            lcqpOASES::OutputStatistics stats;
            lcqp.getOutputStatistics(stats);

            mxArray* iterTotal = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxArray* iterOuter = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxArray* iterSubpr = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxArray* rhoOpt = mxCreateDoubleMatrix(1, 1, mxREAL);
            mxArray* elapsed_time = mxCreateDoubleMatrix(1, 1, mxREAL);

            double* itrTot = mxGetPr(iterTotal);
            double* itrOutr = mxGetPr(iterOuter);
            double* itrSubp = mxGetPr(iterSubpr);
            double* rOpt = mxGetPr(rhoOpt);
            double* elapsed = mxGetPr(elapsed_time);

            itrTot[0] = stats.getIterTotal();
            itrOutr[0] = stats.getIterOuter();
            itrSubp[0] = stats.getSubproblemIter();
            rOpt[0] = stats.getRhoOpt();
            elapsed[0] = elapsed_secs;

            // assign values to struct
            mxSetFieldByNumber(plhs[2], 0, 0, iterTotal);
            mxSetFieldByNumber(plhs[2], 0, 1, iterOuter);
            mxSetFieldByNumber(plhs[2], 0, 2, iterSubpr);
            mxSetFieldByNumber(plhs[2], 0, 3, rhoOpt);
            mxSetFieldByNumber(plhs[2], 0, 4, elapsed_time);
        }
    }

    // Destroy LCQProblem object
    lcqp.~LCQProblem();

    return;
}