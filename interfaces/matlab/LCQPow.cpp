/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2022 by Jonas Hall et al.
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

void readVectors(
    const mxArray** prhs, int nrhs, int nC, double** g,
    double** lbL, double** ubL, double** lbR, double** ubR,
    double** lbA, double** ubA, double** lb, double** ub)
{
    *g = (double*) mxGetPr( prhs[1] );
    *lbL = (double*) mxGetPr( prhs[4] );
    *ubL = (double*) mxGetPr( prhs[5] );
    *lbR = (double*) mxGetPr( prhs[6] );
    *ubR = (double*) mxGetPr( prhs[7] );

    if (nrhs == 10 || (nrhs == 11 && nC == 0)) {
        *lb = (double*) mxGetPr( prhs[8] );
        *ub = (double*) mxGetPr( prhs[9] );
    } else if (nrhs >= 11) {
        *lbA = (double*) mxGetPr( prhs[9] );
        *ubA = (double*) mxGetPr( prhs[10] );

        if (nrhs >= 13) {
            *lb = (double*) mxGetPr( prhs[11] );
            *ub = (double*) mxGetPr( prhs[12] );
        }
    }
}

int LCQPSparse(LCQProblem& lcqp, int nV, int nComp, int nC, int nrhs, const mxArray* prhs[], double* x0, double* y0) {

    if ( !mxIsSparse(prhs[0]) || !mxIsSparse(prhs[2]) || !mxIsSparse(prhs[3]) || (nC > 0 && !mxIsSparse(prhs[8])) )
	{
        mexPrintf("If using the sparse mode, please make sure to provide all matrices in sparse format!\n");
        return 1;
    }

    // Read sparse matrices
    csc* Q = readSparseMatrix(prhs[0], nV, nV);
    csc* L = readSparseMatrix(prhs[2], nComp, nV);
    csc* R = readSparseMatrix(prhs[3], nComp, nV);
    csc* A = NULL;
    if (nC > 0) {
        A = readSparseMatrix(prhs[8], nC, nV);
    }

    // Read vectors
    double* g;
    double* lbL = NULL;
    double* ubL = NULL;
    double* lbR = NULL;
    double* ubR = NULL;
    double* lbA = NULL;
    double* ubA = NULL;
    double* lb = NULL;
    double* ub = NULL;

    // Read all vectors
    readVectors(prhs, nrhs, nC, &g, &lbL, &ubL, &lbR, &ubR, &lbA, &ubA, &lb, &ub);

    // Load data into LCQP object
    LCQPow::ReturnValue ret = lcqp.loadLCQP(
        Q, g, L, R, lbL, ubL, lbR, ubR, A, lbA, ubA, lb, ub, x0, y0
    );

    // Clear sparse specific memory
    LCQPow::Utilities::ClearSparseMat(Q);
    LCQPow::Utilities::ClearSparseMat(L);
    LCQPow::Utilities::ClearSparseMat(R);
    if (A != 0) LCQPow::Utilities::ClearSparseMat(A);

    return ret;
}

int LCQPDense(LCQProblem& lcqp, int nV, int nComp, int nC, int nrhs, const mxArray* prhs[], const double* const x0, const double* const y0)
{
    // Load data
    double* Q = NULL;
    double* g = NULL;
    double* L = NULL;
    double* R = NULL;
    double* lbL = NULL;
    double* ubL = NULL;
    double* lbR = NULL;
    double* ubR = NULL;
    double* lb = NULL;
    double* ub = NULL;
    double* A = NULL;
    double* lbA = NULL;
    double* ubA = NULL;

    // Temporary matrices to keep the column major formats
    double* L_col = NULL;
    double* R_col = NULL;
    double* A_col = NULL;

    Q = (double*) mxGetPr( prhs[0] );
    L_col = (double*) mxGetPr( prhs[2] );
    R_col = (double*) mxGetPr( prhs[3] );

    if (nC > 0) {
        A_col = (double*) mxGetPr( prhs[8] );
    }

    // Read all vectors
    readVectors(prhs, nrhs, nC, &g, &lbL, &ubL, &lbR, &ubR, &lbA, &ubA, &lb, &ub);

    // MATLAB stores in column major format (switch to row major)
    if (L_col != NULL && nComp > 0 && nV > 0) {
        L = new double[nComp*nV];
        colMajorToRowMajor(L_col, L, nComp, nV);
    }
    if (R_col != NULL && nComp > 0 && nV > 0) {
        R = new double[nComp*nV];
        colMajorToRowMajor(R_col, R, nComp, nV);
    }
    if (A_col != NULL && nC > 0 && nV > 0) {
        A = new double[nC*nV];
        colMajorToRowMajor(A_col, A, nC, nV);
    }

    // Load data into LCQP object
    LCQPow::ReturnValue ret = lcqp.loadLCQP(Q, g, L, R, lbL, ubL, lbR, ubR, A, lbA, ubA, lb, ub, x0, y0);

    // Clear A, L, R
    if (A != 0)
        delete[] A;

    if (L != 0)
        delete[] L;

    if (R != 0)
        delete[] R;

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
    int nrhs_min = 8; int nrhs_max = 14;
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
        sprintf(errorMsg, "Invalid input argument: L must be a non-empty double matrix.\n");
        mexErrMsgTxt(errorMsg);
        free(errorMsg);
        return;
    } else {
        nComp = (int)mxGetM(prhs[2]);
    }

    // Get number of linear constraints
    if (nrhs > 11 || ( nrhs == 11 && !mxIsStruct(prhs[10]))) {
        if (mxIsEmpty(prhs[8])) {
            nC = 0;
        } else if (!mxIsDouble(prhs[8])) {
            char *errorMsg = (char*)malloc(100*sizeof(char));
            sprintf(errorMsg, "Invalid input argument: A must be a double matrix.\n");
            mexErrMsgTxt(errorMsg);
            free(errorMsg);
            return;
        } else {
            nC = (int)mxGetM(prhs[8]);
        }
    }

    // Check all dimensions (except for params)
    if (!checkDimensionAndTypeDouble(prhs[0], nV, nV, "Q")) return;
    if (!checkDimensionAndTypeDouble(prhs[1], nV, 1, "g")) return;
    if (!checkDimensionAndTypeDouble(prhs[2], nComp, nV, "L")) return;
    if (!checkDimensionAndTypeDouble(prhs[3], nComp, nV, "R")) return;

    if (nrhs == 10 && !checkDimensionAndTypeDouble(prhs[8], nV, 1, "lb", true)) return;
    if (nrhs == 10 && !checkDimensionAndTypeDouble(prhs[9], nV, 1, "ub", true)) return;

    if (nrhs == 11 && mxIsStruct(prhs[10]) && !checkDimensionAndTypeDouble(prhs[8], nV, 1, "lb", true)) return;
    if (nrhs == 11 && mxIsStruct(prhs[10]) && !checkDimensionAndTypeDouble(prhs[9], nV, 1, "ub", true)) return;

    if (nrhs == 11 && !mxIsStruct(prhs[10]) && !checkDimensionAndTypeDouble(prhs[8], nC, nV, "A", true)) return;
    if (nrhs == 11 && !mxIsStruct(prhs[10]) && !checkDimensionAndTypeDouble(prhs[9], nC, 1, "lbA", true)) return;
    if (nrhs == 11 && !mxIsStruct(prhs[10]) && !checkDimensionAndTypeDouble(prhs[10], nC, 1, "ubA", true)) return;

    if (nrhs >= 12 && !checkDimensionAndTypeDouble(prhs[8], nC, nV, "A", true)) return;
    if (nrhs >= 12 && !checkDimensionAndTypeDouble(prhs[9], nC, 1, "lbA", true)) return;
    if (nrhs >= 12 && !checkDimensionAndTypeDouble(prhs[10], nC, 1, "ubA", true)) return;

    if ((nrhs == 13 || nrhs == 14) && !checkDimensionAndTypeDouble(prhs[11], nV, 1, "lb", true)) return;
    if ((nrhs == 13 || nrhs == 14) && !checkDimensionAndTypeDouble(prhs[12], nV, 1, "ub", true)) return;

    // Check structs
    if (nrhs == 9 && !checkTypeStruct(prhs[8], "params")) return;
    if (nrhs == 11 && !mxIsDouble(prhs[10]) && !checkTypeStruct(prhs[10], "params")) return;
    if (nrhs == 12 && !checkTypeStruct(prhs[11], "params")) return;
    if (nrhs == 14 && !checkTypeStruct(prhs[13], "params")) return;

    // Initial guess variables
    double* x0 = NULL;
    double* y0 = NULL;

    // Create LCQP and options objects
    LCQProblem lcqp((int)nV, (int)nC, (int)nComp);

    // Load settings
    Options options;
    int structIdx = -1;
    if (nrhs == 9 && checkTypeStruct(prhs[8], "params")) { structIdx = 8; }
    if (nrhs == 11 && !mxIsDouble(prhs[10]) && checkTypeStruct(prhs[10], "params"))  { structIdx = 10; }
    if (nrhs == 12 && checkTypeStruct(prhs[11], "params")) { structIdx = 11; }
    if (nrhs == 14 && checkTypeStruct(prhs[13], "params")) { structIdx = 13; }

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
            "maxRho",
            "nDynamicPenalty",
            "etaDynamicPenalty",
            "printLevel",
            "storeSteps",
            "qpSolver"
        };

        mxArray* dual_guess_field;
        bool dual_guess_passed = false;

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

            if ( strcmp(name, "maxRho") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.maxRho")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setMaxRho( (int)fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "nDynamicPenalty") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.nDynamicPenalty")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setNDynamicPenalty( (int)fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "etaDynamicPenalty") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.etaDynamicPenalty")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setEtaDynamicPenalty( fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "printLevel") == 0 ) {
                if (!checkDimensionAndTypeDouble(field, 1, 1, "params.printLevel")) return;

                fld_ptr = (double*) mxGetPr(field);
                options.setPrintLevel( (int)fld_ptr[0] );
                continue;
            }

            if ( strcmp(name, "storeSteps") == 0 ) {
                if (!checkDimensionAndTypeBool(field, 1, 1, "params.storeSteps")) return;

                bool* fld_ptr_bool = (bool*) mxGetPr(field);
                options.setStoreSteps( fld_ptr_bool[0] );
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
                dual_guess_passed = true;
                dual_guess_field = field;
            }
        }

        if (dual_guess_passed) {            
            LCQPow::QPSolver sol = options.getQPSolver();

            if (sol < LCQPow::QPSolver::OSQP_SPARSE) {
                int nDualsIn = nV + nC + 2*nComp;
                if (!checkDimensionAndTypeDouble(dual_guess_field, nDualsIn, 1, "params.y0")) return;
            } else {
                int nDualsIn = nC + 2*nComp;
                 if (!checkDimensionAndTypeDouble(dual_guess_field, nDualsIn, 1, "params.y0")) return;
            }
            
            y0 = (double*) mxGetPr(dual_guess_field);
        }
    }


    // Set options and print them
    lcqp.setOptions( options );

    // For debug sakes
    // printOptions( options );

    // Before running solver, let's allocate solution vectors
    // 1) Primal solution vector
    plhs[0] = mxCreateDoubleMatrix((size_t)nV, 1, mxREAL);
    if (plhs[0] == NULL) {
        mexPrintf("Failed to allocate output of primal solution vector.\n");
        return;
    }
    
    // 2) Dual solution vector (might be too large, but certainly large enough)
    double* yOptTMP = NULL;
    if (nlhs > 1) {
        int maxDualsOut = nC + nV + 2*nComp;
        yOptTMP = new double[maxDualsOut];
    }

    // Point to the output object
    double* xOpt = (double*) mxGetPr(plhs[0]);

    // Start time (time entire c++ code, i.e., result will not include the interface's overhead).
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Sparsity handling
    int ret = 0;
    if (mxIsSparse(prhs[0])) {
        ret = LCQPSparse(lcqp, nV, nComp, nC, nrhs, prhs, x0, y0);
    } else {
        ret = LCQPDense(lcqp, nV, nComp, nC, nrhs, prhs, x0, y0);
    }

    if (ret != 0) {
        mexPrintf("Failed to load LCQP (error code %d).\n", ret);
        return;
    }

    // Run solver
    ret = lcqp.runSolver();

    // Get primal solution 
    lcqp.getPrimalSolution(xOpt);

    // Get dual solution 
    if (nlhs > 1) {
        lcqp.getDualSolution(yOptTMP);
    }

    // Get the statistics
    LCQPow::OutputStatistics stats;
    if (nlhs > 2) {
        lcqp.getOutputStatistics(stats);
    }

    // Stop the timer
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsed_secs = (double)(end - begin).count()/1000.0/1000.0/1000.0;
    
    if (ret != LCQPow::SUCCESSFUL_RETURN) {
        mexPrintf("Failed to solve LCQP (error code: %d).\n", ret);
    }
    
    // Post-Processing
    // If number of duals is actually smaller than nV+nC+2*nComp (e.g. OSQP solver)
    if (nlhs > 1) {

        // Allocate output vector
        int nDualsOut = lcqp.getNumberOfDuals();
        plhs[1] = mxCreateDoubleMatrix((size_t)nDualsOut, 1, mxREAL);

        // Creater pointer to output
        double* yOpt = (double*) mxGetPr(plhs[1]);  

        // Assign dual output
        if (nDualsOut <= 0) {
            mexPrintf("Failed to receive number of dual variables.\n");
            return;
        }
        
        // Copy duals to corrected size
        for (int i = 0; i < nDualsOut; i++) 
            yOpt[i] = yOptTMP[i];
        
        delete[] yOptTMP;

        if (plhs[1] == NULL) {
            mexPrintf("Failed to allocate output of dual solution vector.\n");
            return;
        }
    }

    // Assign the output statistics
    if (nlhs > 2) {
        // assign fieldnames
        int numberStatOutputs = 6;

        if (options.getStoreSteps()) {
            const char* fieldnames[] = {
                "iters_total", "iters_outer", "iters_subproblem", "rho_opt", "elapsed_time", "exit_flag",
                "innerIters", "xSteps", "accumulatedSubproblemIters", "stepLength", "stepSize",
                "statVals", "objVals", "phiVals", "meritVals", "subproblemIters"
            };
            numberStatOutputs += 10;

            // Allocate memory
            plhs[2] = mxCreateStructMatrix(1, 1, numberStatOutputs, fieldnames);
        } else {
            const char* fieldnames[] = {"iters_total", "iters_outer", "iters_subproblem", "rho_opt", "elapsed_time", "exit_flag"};

            // Allocate memory
            plhs[2] = mxCreateStructMatrix(1, 1, numberStatOutputs, fieldnames);
        }

        // assert that allocation went ok
        if (plhs[2] == NULL) {
            mexPrintf("Failed to allocate output of statistics struct.\n");
            return;
        }

        mxArray* iterTotal = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* iterOuter = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* iterSubpr = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* rhoOpt = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* elapsed_time = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxArray* exit_flag = mxCreateDoubleMatrix(1,1, mxREAL);
        mxArray* qp_exit_flag = mxCreateDoubleMatrix(1,1, mxREAL);

        double* itrTot = mxGetPr(iterTotal);
        double* itrOutr = mxGetPr(iterOuter);
        double* itrSubp = mxGetPr(iterSubpr);
        double* rOpt = mxGetPr(rhoOpt);
        double* elapsed = mxGetPr(elapsed_time);
        double* ex_flag = mxGetPr(exit_flag);
        double* qp_ex_flag = mxGetPr(qp_exit_flag);

        itrTot[0] = stats.getIterTotal();
        itrOutr[0] = stats.getIterOuter();
        itrSubp[0] = stats.getSubproblemIter();
        rOpt[0] = stats.getRhoOpt();
        elapsed[0] = elapsed_secs;
        ex_flag[0] = ret;
        qp_ex_flag[0] = stats.getQPSolverExitFlag();

        // assign values to struct
        mxSetFieldByNumber(plhs[2], 0, 0, iterTotal);
        mxSetFieldByNumber(plhs[2], 0, 1, iterOuter);
        mxSetFieldByNumber(plhs[2], 0, 2, iterSubpr);
        mxSetFieldByNumber(plhs[2], 0, 3, rhoOpt);
        mxSetFieldByNumber(plhs[2], 0, 4, elapsed_time);
        mxSetFieldByNumber(plhs[2], 0, 5, exit_flag);
        mxSetFieldByNumber(plhs[2], 0, 6, qp_exit_flag);

        // Tracking values
        if (options.getStoreSteps()) {

            if (stats.getIterTotal() > 0) {
                mxArray* xStepsArr = mxCreateDoubleMatrix((mwSize)(stats.getIterTotal()*nV), 1, mxREAL);
                mxArray* innerItersArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* subproblemItersArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* accuSubproblemItersArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* stepLengthArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* stepSizeArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* statValsArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* objValsArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* phiValsArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);
                mxArray* meritValsArr = mxCreateDoubleMatrix((mwSize)stats.getIterTotal(), 1, mxREAL);

                printf("Iterations: %d\n", stats.getIterTotal());
                printf("Iterations*nV: %d", stats.getIterTotal()*nV);

                double* xSteps = (double*) mxGetPr(xStepsArr);
                double* innerIters = (double*) mxGetPr(innerItersArr);
                double* subproblemIters = (double*) mxGetPr(subproblemItersArr);
                double* accuSubproblemIters = (double*) mxGetPr(accuSubproblemItersArr);
                double* stepLength = (double*) mxGetPr(stepLengthArr);
                double* stepSize = (double*) mxGetPr(stepSizeArr);
                double* statVals = (double*) mxGetPr(statValsArr);
                double* objVals = (double*) mxGetPr(objValsArr);
                double* phiVals = (double*) mxGetPr(phiValsArr);
                double* meritVals = (double*) mxGetPr(meritValsArr);

                std::vector<std::vector<double>> xStepsTMP = stats.getxStepsStdVec( );
                int* innerItersTMP = stats.getInnerIters( );
                int* subproblemItersTMP = stats.getSubproblemIters( );
                int* accuSubproblemItersTMP = stats.getAccuSubproblemIters( );
                double* stepLengthTMP = stats.getStepLength( );
                double* stepSizeTMP = stats.getStepSize( );
                double* statValsTMP = stats.getStatVals( );
                double* objValsTMP = stats.getObjVals( );
                double* phiValsTMP = stats.getPhiVals( );
                double* meritValsTMP = stats.getMeritVals( );

                for (int i = 0; i < stats.getIterTotal(); i++) {
                    for (int j = 0; j < nV; j++) {
                        xSteps[i*nV+j] = xStepsTMP[i][j];
                    }

                    innerIters[i] = innerItersTMP[i];
                    subproblemIters[i] = subproblemItersTMP[i];
                    accuSubproblemIters[i] = accuSubproblemItersTMP[i];
                    stepLength[i] = stepLengthTMP[i];
                    stepSize[i] = stepSizeTMP[i];
                    statVals[i] = statValsTMP[i];
                    objVals[i] = objValsTMP[i];
                    phiVals[i] = phiValsTMP[i];
                    meritVals[i] = meritValsTMP[i];
                }

                mxSetFieldByNumber(plhs[2], 0, 7, xStepsArr);
                mxSetFieldByNumber(plhs[2], 0, 8, subproblemItersArr);
                mxSetFieldByNumber(plhs[2], 0, 9, accuSubproblemItersArr);
                mxSetFieldByNumber(plhs[2], 0, 10, stepLengthArr);
                mxSetFieldByNumber(plhs[2], 0, 11, stepSizeArr);
                mxSetFieldByNumber(plhs[2], 0, 12, statValsArr);
                mxSetFieldByNumber(plhs[2], 0, 13, objValsArr);
                mxSetFieldByNumber(plhs[2], 0, 14, phiValsArr);
                mxSetFieldByNumber(plhs[2], 0, 15, meritValsArr);
                mxSetFieldByNumber(plhs[2], 0, 16, innerItersArr);

                delete[] innerItersTMP;
                delete[] subproblemItersTMP;
                delete[] accuSubproblemItersTMP;
                delete[] stepLengthTMP;
                delete[] stepSizeTMP;
                delete[] statValsTMP;
                delete[] objValsTMP;
                delete[] phiValsTMP;
                delete[] meritValsTMP;
                xStepsTMP.clear();
            }
        }
    }

    return;
}