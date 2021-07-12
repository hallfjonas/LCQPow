/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "LCQProblem.hpp"
using LCQPow::LCQProblem;
using LCQPow::Options;

#include "qpOASES.hpp"

#include <mex.h>
#include <chrono>

bool checkDimensionAndTypeDouble(const mxArray* arr, int m, int n, const char* name, bool allowEmpty = false)
{
    if (allowEmpty && mxIsEmpty(arr)) {
        return true;
    }

    if (!mxIsDouble(arr)) {
        char* errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid type: %s must be of type double.\n", name);
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return false;
    }

    if (mxGetM(arr) != (size_t)m || mxGetN(arr) != (size_t)n) {
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
    mexPrintf(" \n Using LCQPow Options: \n");
    mexPrintf("          rho0: %g \n", options.getInitialPenaltyParameter());
    mexPrintf("          beta: %g \n", options.getPenaltyUpdateFactor());
    mexPrintf("     compl tol: %g \n", options.getComplementarityTolerance());
    mexPrintf("     stati tol: %g \n", options.getStationarityTolerance());
    mexPrintf("      max iter: %d \n", options.getMaxIterations());
    mexPrintf("solve zero pen: %d \n", (int) options.getSolveZeroPenaltyFirst());
    mexPrintf("          prnt: %d \n\n", options.getPrintLevel());
}

csc* readSparseMatrix(const mxArray* mat, int nRow, int nCol)
{
    mwIndex *mat_ir = mxGetIr( mat );
    mwIndex *mat_jc = mxGetJc( mat );
    double *v = (double*)mxGetPr( mat );
    size_t M_nnx = mat_jc[(mwIndex)nCol];
    double* M_data = (double*) malloc(M_nnx*sizeof(double));
    int* M_i = (int*) malloc(M_nnx*sizeof(int));
    int* M_p = (int*) malloc((size_t)(nCol+1)*sizeof(int));
    for (size_t i = 0; i < M_nnx; i++) {
        M_data[i] = v[i];
        M_i[i] = (int) mat_ir[i];
    }

    for (int i = 0; i < nCol+1; i++) {
        M_p[i] = (int) mat_jc[i];
    }

    return LCQPow::Utilities::createCSC(nRow, nCol, M_p[nCol], M_data, M_i, M_p);
}

void readVectors(const mxArray** prhs, int nrhs, int nC, double** g, double** lbA, double** ubA, double** lb, double** ub)
{
    *g = (double*) mxGetPr( prhs[1] );
    if (nrhs == 6 || (nrhs == 7 && nC == 0)) {
        *lb = (double*) mxGetPr( prhs[4] );
        *ub = (double*) mxGetPr( prhs[5] );
    } else if (nrhs >= 7) {
        *lbA = (double*) mxGetPr( prhs[5] );
        *ubA = (double*) mxGetPr( prhs[6] );

        if (nrhs >= 9) {
            *lb = (double*) mxGetPr( prhs[7] );
            *ub = (double*) mxGetPr( prhs[8] );
        }
    }
}

int LCQPSparse(LCQProblem& lcqp, int nV, int nComp, int nC, int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[], double* x0, double* y0) {

    if ( !mxIsSparse(prhs[0]) || !mxIsSparse(prhs[2]) || !mxIsSparse(prhs[3]) || (nC > 0 && !mxIsSparse(prhs[4])) )
	{
        mexPrintf("If using the sparse mode, please make sure to provide all matrices in sparse format!\n");
        return 1;
    }

    // Read sparse matrices
    csc* H = readSparseMatrix(prhs[0], nV, nV);
    csc* S1 = readSparseMatrix(prhs[2], nComp, nV);
    csc* S2 = readSparseMatrix(prhs[3], nComp, nV);
    csc* A = NULL;
    if (nC > 0) {
        A = readSparseMatrix(prhs[4], nC, nV);
    }

    // Read vectors
    double* g;
    double* lbA = NULL;
    double* ubA = NULL;
    double* lb = NULL;
    double* ub = NULL;
    readVectors(prhs, nrhs, nC, &g, &lbA, &ubA, &lb, &ub);

    // Load data into LCQP object
    LCQPow::ReturnValue ret = lcqp.loadLCQP(
        H, g, S1, S2, A, lbA, ubA, lb, ub, x0, y0
    );

    // Clear sparse specific memory
    LCQPow::Utilities::ClearSparseMat(H);
    LCQPow::Utilities::ClearSparseMat(S1);
    LCQPow::Utilities::ClearSparseMat(S2);
    if (A != 0) LCQPow::Utilities::ClearSparseMat(A);

    return ret;
}

int LCQPDense(LCQProblem& lcqp, int nV, int nComp, int nC, int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[], const double* const x0, const double* const y0)
{
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

    // Temporary matrices to keep the column major formats
    double* S1_col = NULL;
    double* S2_col = NULL;
    double* A_col = NULL;

    H = (double*) mxGetPr( prhs[0] );
    S1_col = (double*) mxGetPr( prhs[2] );
    S2_col = (double*) mxGetPr( prhs[3] );

    if (nC > 0) {
        A_col = (double*) mxGetPr( prhs[4] );
    }

    // Read all vectors
    readVectors(prhs, nrhs, nC, &g, &lbA, &ubA, &lb, &ub);

    // MATLAB stores in column major format (switch to row major)
    if (S1_col != NULL && nComp > 0 && nV > 0) {
        S1 = new double[nComp*nV];
        colMajorToRowMajor(S1_col, S1, nComp, nV);
    }
    if (S2_col != NULL && nComp > 0 && nV > 0) {
        S2 = new double[nComp*nV];
        colMajorToRowMajor(S2_col, S2, nComp, nV);
    }
    if (A_col != NULL && nC > 0 && nV > 0) {
        A = new double[nC*nV];
        colMajorToRowMajor(A_col, A, nC, nV);
    }

    // Load data into LCQP object
    LCQPow::ReturnValue ret = lcqp.loadLCQP(H, g, S1, S2, A, lbA, ubA, lb, ub, x0);

    // Clear A, S1, S2
    if (A != 0)
        delete[] A;

    if (S1 != 0)
        delete[] S1;

    if (S2 != 0)
        delete[] S2;

    if (ret != LCQPow::ReturnValue::SUCCESSFUL_RETURN) {
        mexPrintf("Failed to load LCQP.\n");
        return 1;
    }

    return 0;
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

    int nV = 0;
    int nComp = 0;
    int nC = 0;

    // Get number of optimization variables
    if (mxIsEmpty(prhs[0]) || !mxIsDouble(prhs[0])) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid input argument: Hessian must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nV = (int)mxGetM(prhs[0]);
    }

    // Get number of complementarity constraints
    if (mxIsEmpty(prhs[2]) || !mxIsDouble(prhs[2])) {
        char *errorMsg = (char*)malloc(100*sizeof(char));
        sprintf(errorMsg, "Invalid input argument: S1 must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nComp = (int)mxGetM(prhs[2]);
    }

    // Get number of linear constraints
    if (nrhs > 7 || ( nrhs == 7 && !mxIsStruct(prhs[6]))) {
        if (mxIsEmpty(prhs[4])) {
            nC = 0;
        } else if (!mxIsDouble(prhs[4])) {
            char *errorMsg = (char*)malloc(100*sizeof(char));
            sprintf(errorMsg, "Invalid input argument: A must be a double matrix.\n");
            mexErrMsgTxt(errorMsg);
            free(errorMsg);
            return;
        } else {
            nC = (int)mxGetM(prhs[4]);
        }
    }

    // Check all dimensions (except for params)
    if (!checkDimensionAndTypeDouble(prhs[0], nV, nV, "H")) return;
    if (!checkDimensionAndTypeDouble(prhs[1], nV, 1, "g")) return;
    if (!checkDimensionAndTypeDouble(prhs[2], nComp, nV, "S1")) return;
    if (!checkDimensionAndTypeDouble(prhs[3], nComp, nV, "S2")) return;

    if (nrhs == 6 && !checkDimensionAndTypeDouble(prhs[4], nV, 1, "lb", true)) return;
    if (nrhs == 6 && !checkDimensionAndTypeDouble(prhs[5], nV, 1, "ub", true)) return;

    if (nrhs == 7 && mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[4], nV, 1, "lb", true)) return;
    if (nrhs == 7 && mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[5], nV, 1, "ub", true)) return;

    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[4], nC, nV, "A", true)) return;
    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[5], nC, 1, "lbA", true)) return;
    if (nrhs == 7 && !mxIsStruct(prhs[6]) && !checkDimensionAndTypeDouble(prhs[6], nC, 1, "ubA", true)) return;

    if (nrhs >= 8 && !checkDimensionAndTypeDouble(prhs[4], nC, nV, "A", true)) return;
    if (nrhs >= 8 && !checkDimensionAndTypeDouble(prhs[5], nC, 1, "lbA", true)) return;
    if (nrhs >= 8 && !checkDimensionAndTypeDouble(prhs[6], nC, 1, "ubA", true)) return;

    if ((nrhs == 9 || nrhs == 10) && !checkDimensionAndTypeDouble(prhs[7], nV, 1, "lb", true)) return;
    if ((nrhs == 9 || nrhs == 10) && !checkDimensionAndTypeDouble(prhs[8], nV, 1, "ub", true)) return;

    // Check structs
    if (nrhs == 5 && !checkTypeStruct(prhs[4], "params")) return;
    if (nrhs == 7 && !mxIsDouble(prhs[6]) && !checkTypeStruct(prhs[6], "params")) return;
    if (nrhs == 8 && !checkTypeStruct(prhs[7], "params")) return;
    if (nrhs == 10 && !checkTypeStruct(prhs[9], "params")) return;

    // Initial guess variables
    double* x0 = NULL;
    double* y0 = NULL;

    // Create LCQP and options objects
    LCQProblem lcqp((int)nV, (int)nC, (int)nComp);

    // Load settings
    Options options;
    int structIdx = -1;
    if (nrhs == 5 && checkTypeStruct(prhs[4], "params")) { structIdx = 4; }
    if (nrhs == 7 && !mxIsDouble(prhs[6]) && checkTypeStruct(prhs[6], "params"))  { structIdx = 6; }
    if (nrhs == 8 && checkTypeStruct(prhs[7], "params")) { structIdx = 7; }
    if (nrhs == 10 && checkTypeStruct(prhs[9], "params")) { structIdx = 9; }

    if (structIdx != -1) {

        const char* params_fieldnames[] = {
            "x0",
            "y0",
            "stationarityTolerance",
            "complementarityTolerance",
            "initialPenaltyParameter",
            "penaltyUpdateFactor",
            "solveZeroPenaltyFirst",
            "maxIterations",
            "printLevel",
            "qpSolver"
        };

        for (auto name : params_fieldnames) {
            mxArray* field = mxGetField(prhs[structIdx], 0, name);
            double* fld_ptr;

            if (field == NULL)
                continue;

            if ( strcmp(name, "stationarityTolerance") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.stationarityTolerance")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setStationarityTolerance( fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "complementarityTolerance") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.complementarityTolerance")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setComplementarityTolerance( fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "initialPenaltyParameter") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.initialPenaltyParameter")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setInitialPenaltyParameter( fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "penaltyUpdateFactor") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.penaltyUpdateFactor")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setPenaltyUpdateFactor( fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "solveZeroPenaltyFirst") == 0 ) {
                if (!checkDimensionAndTypeBool(field, 1, 1, "params.solveZeroPenaltyFirst")) return;

                bool* fld_ptr_bool = (bool*) mxGetPr(field);
                options.setSolveZeroPenaltyFirst( fld_ptr_bool[0] );
                continue;
            }

            if ( strcmp(name, "maxIterations") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.maxIterations")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setMaxIterations( (int)fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "printLevel") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.printLevel")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setPrintLevel( (int)fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "qpSolver") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.qpSolver")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setQPSolver( (int)fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "x0") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, nV, 1, "params.x0")) return;

                x0 = (double*) mxGetPr(field);
                continue;
            }

            if ( strcmp(name, "y0") == 0 ) {
                // TODO: This is just a heuristic right now (OSQP does not take box constraints, but we don't know th QP solver yet...)
                int nDualsIn1 = nV + nC + 2*nComp;
                int nDualsIn2 = nC + 2*nComp;

                if (!checkDimensionAndTypeDouble(field, nDualsIn1, 1, "params.y0") && !checkDimensionAndTypeDouble(field, nDualsIn2, 1, "params.y0")) return;

                y0 = (double*) mxGetPr(field);
                continue;
            }
        }
    }

    // Set options and print them
    lcqp.setOptions( options );

    // For debug sakes
    // printOptions( options );

    // Sparsity checks
    int ret = 0;
    if (mxIsSparse(prhs[0]) || mxIsSparse(prhs[2]) || mxIsSparse(prhs[3])|| (nC > 0 && mxIsSparse(prhs[4]))) {
        ret = LCQPSparse(lcqp, nV, nComp, nC, nlhs, plhs, nrhs, prhs, x0, y0);
    } else {
        ret = LCQPDense(lcqp, nV, nComp, nC, nlhs, plhs, nrhs, prhs, x0, y0);
    }

    if (ret != 0) {
        mexPrintf("Failed to load LCQP (error code %d).\n", ret);
        return;
    }

    // Start time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Run solver
    ret = lcqp.runSolver();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsed_secs = (double)(end - begin).count()/1000.0/1000.0/1000.0;

    if (ret != LCQPow::SUCCESSFUL_RETURN) {
        mexPrintf("Failed to solve LCQP (error code: %d).\n", ret);
    }

    // Succeeded to solve LCQP. Obtaining solution vector
    // 1) Primal solution vector
    plhs[0] = mxCreateDoubleMatrix((size_t)nV, 1, mxREAL);
    if (plhs[0] == NULL) {
        mexPrintf("Failed to allocate output of primal solution vector.\n");
        return;
    }

    // Point to the output object
    double* xOpt = (double*) mxGetPr(plhs[0]);
    lcqp.getPrimalSolution(xOpt);

    // 2) Dual solution vector
    if (nlhs > 1) {
        int nDualsOut = lcqp.getNumerOfDuals();

        if (nDualsOut <= 0) {
            mexPrintf("Failed to receive number of dual variables.\n");
            return;
        }

        plhs[1] = mxCreateDoubleMatrix((size_t)nDualsOut, 1, mxREAL);
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
        int numberStatOutputs = 6;

        const char* fieldnames[] = {"iters_total", "iters_outer", "iters_subproblem", "rho_opt", "elapsed_time", "exit_flag"};

        // Allocate memory
        plhs[2] = mxCreateStructMatrix(1,1,numberStatOutputs, fieldnames);

        // assert that allocation went ok
        if (plhs[2] == NULL) {
            mexPrintf("Failed to allocate output of statistics struct.\n");
            return;
        }

        // Get the statistics
        LCQPow::OutputStatistics stats;
        lcqp.getOutputStatistics(stats);

        mxArray* iterTotal = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* iterOuter = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* iterSubpr = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* rhoOpt = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* elapsed_time = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* exit_flag = mxCreateDoubleMatrix(1,1, mxREAL);

        double* itrTot = mxGetPr(iterTotal);
        double* itrOutr = mxGetPr(iterOuter);
        double* itrSubp = mxGetPr(iterSubpr);
        double* rOpt = mxGetPr(rhoOpt);
        double* elapsed = mxGetPr(elapsed_time);
        double* ex_flag = mxGetPr(exit_flag);

        itrTot[0] = stats.getIterTotal();
        itrOutr[0] = stats.getIterOuter();
        itrSubp[0] = stats.getSubproblemIter();
        rOpt[0] = stats.getRhoOpt();
        elapsed[0] = elapsed_secs;
        ex_flag[0] = ret;

        // assign values to struct
        mxSetFieldByNumber(plhs[2], 0, 0, iterTotal);
        mxSetFieldByNumber(plhs[2], 0, 1, iterOuter);
        mxSetFieldByNumber(plhs[2], 0, 2, iterSubpr);
        mxSetFieldByNumber(plhs[2], 0, 3, rhoOpt);
        mxSetFieldByNumber(plhs[2], 0, 4, elapsed_time);
        mxSetFieldByNumber(plhs[2], 0, 5, exit_flag);
    }

    return;
}