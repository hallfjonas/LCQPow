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


#include "Utilities.hpp"
#include <iostream>
#include <vector>
#include <osqp.h>

namespace lcqpOASES {


    Options::Options( )
    {
        setToDefault( );
    }


    Options::Options( const Options& rhs)
    {
        copy( rhs );
    }


    Options::~Options( )
    { }


    Options& Options::operator=( const Options& rhs )
    {
        if ( this != &rhs )
        {
            copy( rhs );
        }

        return *this;
    }


    void Options::copy( const Options& rhs ) {
        stationarityTolerance = rhs.stationarityTolerance;
        complementarityTolerance = rhs.complementarityTolerance;
        initialPenaltyParameter = rhs.initialPenaltyParameter;
        penaltyUpdateFactor = rhs.penaltyUpdateFactor;
        solveZeroPenaltyFirst = rhs.solveZeroPenaltyFirst;
        maxIterations = rhs.maxIterations;
        printLevel = rhs.printLevel;
    }


    double Options::getStationarityTolerance( ) {
        return stationarityTolerance;
    }


    void Options::setStationarityTolerance( double val ) {
        if (val <= Utilities::EPS)
            throw INVALID_STATIONARITY_TOLERANCE;

        stationarityTolerance = val;
    }


    double Options::getComplementarityTolerance( ) {
        return complementarityTolerance;
    }


    void Options::setComplementarityTolerance( double val ) {
        if (val <= Utilities::EPS)
            throw INVALID_COMPLEMENTARITY_TOLERANCE;

        complementarityTolerance = val;
    }


    double Options::getInitialPenaltyParameter( ) {
        return initialPenaltyParameter;
    }


    void Options::setInitialPenaltyParameter( double val ) {
        if (val <= Utilities::ZERO)
            throw INVALID_INITIAL_PENALTY_VALUE;

        initialPenaltyParameter = val;
    }


    double Options::getPenaltyUpdateFactor( ) {
        return penaltyUpdateFactor;
    }


    void Options::setPenaltyUpdateFactor( double val ) {
        if (val <= 1)
            throw INVALID_PENALTY_UPDATE_VALUE;

        penaltyUpdateFactor = val;
    }


    bool Options::getSolveZeroPenaltyFirst( ) {
        return solveZeroPenaltyFirst;
    }


    void Options::setSolveZeroPenaltyFirst( bool val ) {
        solveZeroPenaltyFirst = val;
    }


    int Options::getMaxIterations( ) {
        return maxIterations;
    }


    void Options::setMaxIterations( int val ) {
        if (val <= 0)
            throw INVALID_MAX_ITERATIONS_VALUE;

        maxIterations = val;
    }


    PrintLevel Options::getPrintLevel( ) {
        return printLevel;
    }


    void Options::setPrintLevel( PrintLevel val ) {
        printLevel = val;
    }


    void Options::setPrintLevel( int val ) {

        if (val < PrintLevel::NONE || val > PrintLevel::SUBPROBLEM_SOLVER_ITERATES)
            throw INVALID_PRINT_LEVEL_VALUE;

        printLevel = (PrintLevel)val;
    }


    void Options::setToDefault( ) {
        complementarityTolerance = 1.0e3 * Utilities::EPS;
        stationarityTolerance  = 1.0e3 * Utilities::EPS*1000;
        initialPenaltyParameter = 0.01;
    	penaltyUpdateFactor  = 2.0;

        solveZeroPenaltyFirst = true;

        maxIterations = 1000;

        printLevel = PrintLevel::INNER_LOOP_ITERATES;
    }


    void Utilities::MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                C[i*p + j] = 0;
                for (int k = 0; k < n; k++) {
                    C[i*p + j] += A[i*n + k]*B[k*p + j];
                }
            }
        }
    }


    void Utilities::TransponsedMatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p) {

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                C[i*p + j] = 0;

                for (int k = 0; k < m; k++) {
                    C[i*p+j] += A[k*n + i]*B[k*p + j];
                }
            }
        }
    }


    void Utilities::MatrixSymmetrizationProduct(const double* const A, const double* const B, double* C, int m, int n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                C[i*n + j] = 0;
                for (int k = 0; k < m; k++) {
                    C[i*n + j] += A[k*n + i]*B[k*n + j] + B[k*n + i]*A[k*n + j];
                }

                // Make symmetric
                C[j*n + i] = C[i*n + j];
            }
        }
    }


    void Utilities::AffineLinearTransformation(const double alpha, const double* const A, const double* const b, const double* const c, double* d, int m, int n) {
        for (int i = 0; i < m; i++) {

            double tmp = 0;
            for (int k = 0; k < n; k++) {
                tmp += A[i*n + k]*b[k];
            }

            d[i] = alpha*tmp + c[i];
        }
    }


    void Utilities::WeightedMatrixAdd(const double alpha, const double* const A, const double beta, const double* const B, double* C, int m, int n) {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i*n + j] = alpha*A[i*n+j] + beta*B[i*n+j];
    }


    void Utilities::WeightedVectorAdd(const double alpha, const double* const a, const double beta, const double* const b, double* c, int m) {
        WeightedMatrixAdd(alpha, a, beta, b, c, m, 1);
    }


    double Utilities::QuadraticFormProduct(const double* const Q, const double* const p, int m) {
        double ret = 0;
        for (int i = 0; i < m; i++) {
            double tmp = 0;
            for (int j = 0; j < m; j++)
                tmp += Q[i*m + j]*p[j];

            ret += tmp*p[i];
        }

        return ret;
    }


    double Utilities::DotProduct(const double* const a, const double* const b, int m) {
        double ret = 0;
        for (int i = 0; i < m; i++)
            ret += a[i]*b[i];

        return ret;
    }


    double Utilities::MaxAbs(const double* const a, int m) {
        double mx = 0;
        double mn = 0;

        for (int i = 0; i < m; i++) {
            if (a[i] > mx)
                mx = a[i];
            else if (a[i] < mn)
                mn = a[i];
        }

        return Utilities::getMax(mx, -mn);
    }


    ReturnValue Utilities::readFromFile( int* data, int n, const char* datafilename )
    {
        int i;
        FILE* datafile;

        /* 1) Open file. */
        if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
        {
            fclose( datafile );
            return UNABLE_TO_READ_FILE;
        }

        /* 2) Read data from file. */
        for( i=0; i<n; ++i )
        {
            if ( fscanf( datafile, "%d\n", &(data[i]) ) == 0 )
            {
                fclose( datafile );
                return UNABLE_TO_READ_FILE;
            }
        }

        /* 3) Close file. */
        fclose( datafile );

        return SUCCESSFUL_RETURN;
    }


    ReturnValue Utilities::readFromFile( double* data, int n, const char* datafilename )
    {
        int i;
        FILE* datafile;

        /* 1) Open file. */
        if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
        {
            fclose( datafile );
            return UNABLE_TO_READ_FILE;
        }

        /* 2) Read data from file. */
        for( i=0; i<n; ++i )
        {
            if ( fscanf( datafile, "%lf\n", &(data[i]) ) == 0 )
            {
                fclose( datafile );
                return UNABLE_TO_READ_FILE;
            }
        }

        /* 3) Close file. */
        fclose( datafile );

        return SUCCESSFUL_RETURN;
    }


    ReturnValue Utilities::writeToFile( double* data, int n, const char* datafilename )
    {
        int i;
        FILE* datafile;

        /* 1) Open file. */
        if ( ( datafile = fopen( datafilename, "w" ) ) == 0 )
        {
            fclose( datafile );
            return UNABLE_TO_READ_FILE;
        }

        /* 2) Read data from file. */
        for( i=0; i<n; ++i )
        {
            if ( fprintf( datafile, "%f\n", data[i] ) == 0 )
            {
                fclose( datafile );
                return UNABLE_TO_READ_FILE;
            }
        }

        /* 3) Close file. */
        fclose( datafile );

        return SUCCESSFUL_RETURN;
    }


    void Utilities::printMatrix(const double* const A, int m, int n, const char* const name)
    {
        printf("Printing matrix %s:\n", name);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                printf("%.5f ", A[i*n + j]);


            printf("\n");
        }

        printf("\n");
    }


    void Utilities::printMatrix(const int* const A, int m, int n, const char* const name)
    {
        printf("Printing matrix %s:\n", name);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                printf("%d ", A[i*n + j]);


            printf("\n");
        }

        printf("\n");
        fflush(stdout);
    }


    void Utilities::printMatrix(const csc* A, const char* const name)
    {
        if (A->m*A->n <= 0)
            return;

        // Get dense representation
        double* dense = new double[A->m*A->n]();
        Utilities::csc_to_dns(A, dense, A->m, A->n);

        // Print the dense matrix
        Utilities::printMatrix(dense, A->m, A->n, name);

        // Clear memory
        delete[] dense;
    }


    void Utilities::printStep(double* xk, double* pk, double* xk_new, double alpha, int nV)
    {
        printf("Printing Step:\n");

        for (int i = 0; i < nV; i++)
            printf("%.2f + %.2f * %.2f = %.2f \n", xk[i], alpha, pk[i], xk_new[i]);

        printf("\n");
    }


    void Utilities::printBounds(double* lb, double* xk, double* ub, int m)
    {
        printf("Printing box constraints:\n");

        for (int i = 0; i < m; i++)
            printf("%.2f <= %.2f <= %.2f \n", lb[i], xk[i], ub[i]);

        printf("\n");
    }


    ReturnValue Utilities::csc_to_dns(const csc* const sparse, double* full, int m, int n)
    {
        for (int j = 0; j < n; j++) {
			for (int i = sparse->p[j]; i < sparse->p[j+1]; i++) {
                // Reached final element
                if (i == sparse->nzmax) {
                    return ReturnValue::SUCCESSFUL_RETURN;
                }

                // Ensure validity of index
                if (sparse->i[i]*n + j >= m*n || sparse->i[i]*n + j < 0) {
                    return MessageHandler::PrintMessage( ReturnValue::INDEX_OUT_OF_BOUNDS );
                }

				full[sparse->i[i]*n + j] = sparse->x[i];
			}
		}

        return ReturnValue::SUCCESSFUL_RETURN;
    }


    csc* Utilities::dns_to_csc(const double* const full, int m, int n)
    {
        csc* sparse = (csc*) c_malloc(sizeof(csc));

        std::vector<double> H_data;
        std::vector<int> H_i;
        sparse->m = m;
        sparse->n = n;
        sparse->nz = -1;
        sparse->p = new int[n+1]();

        for (int i = 0; i < n; i++) {
            // Begin column pointer with previous value
            sparse->p[i+1] = sparse->p[i];

            for (int j = 0; j < m; j++) {
                if (full[j*n + i] > 0 || full[j*n + i] < 0) {
                    H_data.push_back(full[j*n + i]);
                    H_i.push_back(j);
                    sparse->p[i+1]++;
                }
            }
        }

        // Final pointer should point to equal of elemtns
        int nnx = (int) H_i.size();
        sparse->p[n] = nnx;

        sparse->nzmax = nnx;
        sparse->i = new int[nnx];
        sparse->x = new double[nnx];

        for (int i = 0; i < nnx; i++) {
            sparse->i[i] = H_i[(size_t)i];
            sparse->x[i] = H_data[(size_t)i];
        }

        return sparse;
    }


    double Utilities::getAbs(double x)
    {
        #ifdef __NO_FMATH__
        return (x>=0.0) ? x : -x;
        #else
        return fabs(x);
        #endif
    }


    bool Utilities::isEqual(double x, double y, double TOL)
    {
        if ( getAbs(x-y) <= TOL )
            return true;
        else
            return false;
    }


    bool Utilities::isEqual(double x, double y)
    {
        return isEqual(x, y, Utilities::ZERO);
    }


    bool Utilities::isZero(double x, double TOL)
    {
        if ( getAbs(x) <= TOL )
            return true;
        else
            return false;
    }


    bool Utilities::isZero(double x)
    {
        return isZero(x, Utilities::ZERO);
    }


    double Utilities::getSign(double arg)
    {
        if ( arg >= 0.0 )
            return 1.0;
        else
            return -1.0;
    }


    int Utilities::getMax(int x, int y)
    {
        return (y<x) ? x : y;
    }


    int Utilities::getMin(int x, int y)
    {
        return (y>x) ? x : y;
    }


    double Utilities::getMax(double x, double y)
    {
        return (y<x) ? x : y;
    }


    double Utilities::getMin(double x, double y)
    {
        return (y>x) ? x : y;
    }


    OutputStatistics::OutputStatistics( ) { }


    OutputStatistics& OutputStatistics::operator=( const OutputStatistics& rhs )
    {
        iter_total = rhs.iter_total;
        iter_outer = rhs.iter_outer;
        subproblem_iter = rhs.subproblem_iter;
        rho_opt = rhs.rho_opt;
        status = rhs.status;

        return *this;
    }


    void OutputStatistics::reset( )
    {
        iter_total = 0;
        iter_outer = 0;
        subproblem_iter = 0;
        rho_opt = 0.0;
        status = PROBLEM_NOT_SOLVED;
    }


    ReturnValue OutputStatistics::updateIterTotal( int delta_iter )
    {
        if (delta_iter < 0) return INVALID_TOTAL_ITER_COUNT;

        iter_total += delta_iter;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateIterOuter( int delta_iter )
    {
        if (delta_iter < 0) return INVALID_TOTAL_OUTER_ITER;

        iter_outer += delta_iter;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateSubproblemIter( int delta_iter )
    {
        if (delta_iter < 0) return IVALID_SUBPROBLEM_ITER;

        subproblem_iter += delta_iter;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateRhoOpt( double _rho )
    {
        if (_rho <= 0) return INVALID_RHO_OPT;

        rho_opt = _rho;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateSolutionStatus( AlgorithmStatus _status )
    {
        status = _status;
        return SUCCESSFUL_RETURN;
    }


    int OutputStatistics::getIterTotal( ) const
    {
        return iter_total;
    }


    int OutputStatistics::getIterOuter( ) const
    {
        return iter_outer;
    }


    int OutputStatistics::getSubproblemIter( ) const
    {
        return subproblem_iter;
    }


    double OutputStatistics::getRhoOpt( ) const
    {
        return rho_opt;
    }


    AlgorithmStatus OutputStatistics::getSolutionStatus( ) const
    {
        return status;
    }


    ReturnValue MessageHandler::PrintMessage( ReturnValue ret) {

        switch (ret) {
            case SUCCESSFUL_RETURN:
                break;

            case NOT_YET_IMPLEMENTED:
                printf("This method has not yet been implemented.\n");
                break;

            case LCQPOBJECT_NOT_SETUP:
                printf("ERROR: The LCQP object has not been set up correctly.\n");
                break;

            case INDEX_OUT_OF_BOUNDS:
                printf("ERROR: Index out of bounds.\n");
                break;

            case SUBPROBLEM_SOLVER_ERROR:
                printf("ERROR: The subproblem solver produced an error.\n");
                break;

            case UNABLE_TO_READ_FILE:
                printf("ERROR: Unable to read file.\n");
                break;

            case MAX_ITERATIONS_REACHED:
                printf("ERROR: Maximum number of iterations reached.\n");
                break;

            case INITIAL_SUBPROBLEM_FAILED:
                printf("ERROR: Failed to solve initial QP.\n");
                break;

            case INVALID_ARGUMENT:
                printf("ERROR: Invalid argument passed.\n");
                break;

            case INVALID_NUMBER_OF_OPTIM_VARS:
                printf("ERROR: Invalid optimization variable dimension passed (required to be > 0).\n");
                break;

            case INVALID_NUMBER_OF_COMP_VARS:
                printf("ERROR: Invalid complementarity dimension passed (required to be > 0).\n");
                break;

            case INVALID_NUMBER_OF_CONSTRAINT_VARS:
                printf("ERROR: Invalid number of optimization variables passed (required to be >= 0).\n");
                break;

            case INVALID_COMPLEMENTARITY_TOLERANCE:
                printf("ERROR: Invalid argument passed (complementarity tolerance).\n");
                break;

            case INVALID_INITIAL_PENALTY_VALUE:
                printf("ERROR: Invalid argument passed (initial penalty value).\n");
                break;

            case INVALID_PENALTY_UPDATE_VALUE:
                printf("ERROR: Invalid argument passed (penalty update value).\n");
                break;

            case INVALID_MAX_ITERATIONS_VALUE:
                printf("ERROR: Invalid argument passed (maximum iterations).\n");
                break;

            case INVALID_STATIONARITY_TOLERANCE:
                printf("ERROR: Invalid argument passed (stationarity tolerance).\n");
                break;

            case INVALID_RELAX_OPTIONS_TOLERANCE:
                printf("ERROR: Invalid argument passed (relax optiopns tolerance).\n");
                break;

            case INVALID_INDEX_POINTER:
                printf("ERROR: Invalid index pointer passed in csc format.\n");
                break;

            case INVALID_INDEX_ARRAY:
                printf("ERROR: Invalid index array passed in csc format.\n");
                break;

            case INVALID_OSQP_BOX_CONSTRAINTS:
                printf("ERROR: Invalid constraints passed to OSQP solver: This solver does not handle box constraints, please pass them through linear constraints.\n");
                break;

            case INVALID_TOTAL_ITER_COUNT:
                printf("ERROR: Invalid total number of iterations delta passed to output statistics (must be non-negative integer).\n");
                break;

            case INVALID_TOTAL_OUTER_ITER:
                printf("ERROR: Invalid total number of outer iterations delta passed to output statistics (must be non-negative integer).\n");
                break;

            case IVALID_SUBPROBLEM_ITER:
                printf("ERROR: Invalid total number of subproblem solver iterates delta passed to output statistics (must be non-negative integer).\n");
                break;

            case INVALID_RHO_OPT:
                printf("ERROR: Invalid rho value at solution passed to output statistics (must be positive double).\n");
                break;

            case INVALID_PRINT_LEVEL_VALUE:
                printf("ERROR: Invalid integer to be parsed to print level passed (must be in range of enum).\n");
                break;
        }

        fflush(stdout);

        return ret;
    }


    AlgorithmStatus MessageHandler::PrintSolution( AlgorithmStatus algoStat ) {

        if ( algoStat != PROBLEM_NOT_SOLVED)
            printf("\n\n#################################\n");

        switch (algoStat) {
            case PROBLEM_NOT_SOLVED:
                printf("The LCQP has not been solved.\n");
                break;

            case W_STATIONARY_SOLUTION:
                printf("## W-Stationary solution found ##\n");
                break;

            case C_STATIONARY_SOLUTION:
                printf("## C-Stationary solution found ##\n");
                break;

            case M_STATIONARY_SOLUTION:
                printf("## M-Stationary solution found ##\n");
                break;

            case S_STATIONARY_SOLUTION:
                printf("## S-Stationary solution found ##\n");
                break;
        }

        if ( algoStat != PROBLEM_NOT_SOLVED)
            printf("#################################\n\n");

        return algoStat;
    }
}

