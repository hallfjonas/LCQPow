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

namespace lcqpOASES {
    /*
    *	O p t i o n s
    */
    Options::Options( )
    {
        setToDefault( );
    }


    /*
    *	O p t i o n s
    */
    Options::Options( const Options& rhs)
    {
        copy( rhs );
    }


    /*
    *	~ O p t i o n s
    */
    Options::~Options( )
    { }


    /*
    *	o p e r a t o r =
    */
    Options& Options::operator=( const Options& rhs )
    {
        if ( this != &rhs )
        {
            copy( rhs );
        }

        return *this;
    }


    /*
     *   c o p y
     */
    void Options::copy( const Options& rhs ) {
        stationarityTolerance = rhs.stationarityTolerance;
        complementarityTolerance = rhs.complementarityTolerance;
        initialComplementarityPenalty = rhs.initialComplementarityPenalty;
        complementarityPenaltyUpdate = rhs.complementarityPenaltyUpdate;
        solveZeroPenaltyFirst = rhs.solveZeroPenaltyFirst;

        maxOuterIterations = rhs.maxOuterIterations;
        maxInnerIterations = rhs.maxInnerIterations;

        printLvl = rhs.printLvl;
        qpSubSolver = rhs.qpSubSolver;
    }


    /*
     *  e n s u r e C o n s i s t e n c y
     */
    returnValue Options::ensureConsistency( ) {
        
        if (complementarityPenaltyUpdate <= 1)
            throw INVALID_PENALTY_UPDATE_VALUE;

        if (complementarityTolerance < Utilities::EPS)
            throw INVALID_COMPLEMENTARITY_TOLERANCE;

        if (initialComplementarityPenalty <= 0)
            throw INVALID_INITIAL_PENALTY_VALUE;

        if (maxOuterIterations <= 0)
            throw INVALID_MAX_OUTER_ITERATIONS_VALUE;

        if (maxInnerIterations <= 0)
            throw INVALID_MAX_INNER_ITERATIONS_VALUE;

        return SUCCESSFUL_RETURN;
    }


    /*
     *   s e t T o D e f a u l t
     */
    void Options::setToDefault( ) {
        complementarityTolerance = 1.0e3 * Utilities::EPS;
        stationarityTolerance  = 1.0e3 * Utilities::EPS;
        initialComplementarityPenalty = 0.01;
    	complementarityPenaltyUpdate  = 2.0;

        solveZeroPenaltyFirst = true;

        maxOuterIterations = 100;
        maxInnerIterations = 1000;



        qpSubSolver = qpSubproblemSolver::QPOASES;

        printLvl = printLevel::INNER_LOOP_ITERATES;

    }

    /*
     *   M a t r i x M u l t i p l i c a t i o n
     */
    void Utilities::MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                C[i*m + j] = 0;
                for (int k = 0; k < n; k++) {
                    C[i*m + j] += A[i*n + k]*B[k*p + j]; 
                }                
            }
        }
    }


    /*
     *  M a t r i x S y m m e t r i z a t i o n P r o d u c t
     */
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

    /*
     *   A f f i n e L i n e a r T r a n s f o r m a t i o n
     */
    void Utilities::AffineLinearTransformation(const double alpha, const double* const A, const double* const b, const double* const c, double* d, int m, int n) {
        for (int i = 0; i < m; i++) {
            
            double tmp = 0;
            for (int k = 0; k < n; k++) {
                tmp += A[i*n + k]*b[k]; 
            }                

            d[i] = alpha*tmp + c[i];
        }
    }

    /*
     *	r e a d F r o m F i l e
     */
    void Utilities::WeightedMatrixAdd(const double alpha, const double* const A, const double beta, const double* const B, double* C, int m, int n) {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i*n + j] = alpha*A[i*n+j] + beta*B[i*n+j];
    }

    /*
     *	r e a d F r o m F i l e
     */
    void Utilities::WeightedVectorAdd(const double alpha, const double* const a, const double beta, const double* const b, double* c, int m) {
        WeightedMatrixAdd(alpha, a, beta, b, c, m, 1);
    }
        

    /*
     *	Q u a d r a t i c F o r m P r o d u c t
     */
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

    /*
     *	Q u a d r a t i c F o r m P r o d u c t
     */
    double Utilities::DotProduct(const double* const a, const double* const b, int m) {
        double ret = 0;
        for (int i = 0; i < m; i++)
            ret += a[i]*b[i];

        return ret;
    }

    /*
     *	Q u a d r a t i c F o r m P r o d u c t
     */
    double Utilities::MaxAbs(const double* const a, int m) {
        double max = 0;
        double min = 0;

        for (int i = 0; i < m; i++) {
            if (a[i] > max)
                max = a[i];
            else if (a[i] < min)
                min = a[i];            
        }

        return std::max(max, -min);
    }

    /*
     *	r e a d F r o m F i l e
     */
    returnValue Utilities::readFromFile( int* data, int n, const char* datafilename )
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


    /*
     *	r e a d F r o m F i l e
     */
    returnValue Utilities::readFromFile( double* data, int n, const char* datafilename )
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
}

