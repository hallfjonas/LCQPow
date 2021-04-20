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


#ifndef LCQPOASES_UTILITIES_HPP
#define LCQPOASES_UTILITIES_HPP

#include "Types.hpp"

#include <osqp.h>

namespace lcqpOASES {

    enum returnValue {
        // Special values
        NOT_YET_IMPLEMENTED = -1,                       /**< Not yet implemented (internal use only). */
        SUCCESSFUL_RETURN = 0,						    /**< Successful return. */

        // Invalid arguments
        INVALID_ARGUMENT = 100,                         /**< Invalid argument. */
        INVALID_PENALTY_UPDATE_VALUE = 101,             /**< Invalid penalty update value. Needs to be > 1. */
        INVALID_COMPLEMENTARITY_TOLERANCE = 102,        /**< Invalid complementarity tolerance. Must be no smaller than machine precision. */
        INVALID_INITIAL_PENALTY_VALUE = 103,            /**< Invalid initial penalty parameter. Must be positive. */
        INVALID_MAX_OUTER_ITERATIONS_VALUE = 104,       /**< Invalid number of maximal outer iterations. Must be a positive integer */
        INVALID_MAX_INNER_ITERATIONS_VALUE = 105,       /**< Invalid number of maximal inner iterations. Must be a positive integer. */
        INVALID_NUMBER_OF_OPTIM_VARS = 106,             /**< Invalid number of optimization variables. Must be a positive integer. */
        INVALID_NUMBER_OF_COMP_VARS = 107,              /**< Invalid number of complementarity constraints. Must be a positive integer. */
        INVALID_NUMBER_OF_CONSTRAINT_VARS = 108,        /**< Invalid number of linear constraints. Must be a non-negative integer. */
        INVALID_RELAX_OPTIONS_TOLERANCE = 109,          /**< Invalid number of active set changes to switch to precision mode. Must be a positive integer. */

        // Algorithmic errors
        MAX_OUTER_ITERATIONS_REACHED = 200,             /**< Maximum number of outer iterations reached. */
        MAX_INNER_ITERATIONS_REACHED = 201,             /**< Maximum number of inner iterations reached. */
        INITIAL_SUBPROBLEM_FAILED = 202,                /**< Failed to solve the initial QP. */
        SUBPROBLEM_SOLVER_ERROR = 203,                  /**< An error occured in the subproblem solver. */

        // Generic errors
        LCQPOBJECT_NOT_SETUP = 300,                     /**< Constructor has not been called. */
        INDEX_OUT_OF_BOUNDS = 301,                      /**< Index out of bounds. */
        UNABLE_TO_READ_FILE = 302,                      /**< Unable to read a file. */

        // Sparse matrices
        INVALID_INDEX_POINTER = 400,                    /**< Invalid index pointer for a csc matrix. */
        INVALID_INDEX_ARRAY = 401                       /**< Invalid index array for a csc matrix. */
    };

    enum algorithmStatus {
        PROBLEM_NOT_SOLVED = 0,                         /**< The problem was not solved. */
        W_STATIONARY_SOLUTION = 1,                      /**< The solution corresponds to a weakly stationary point. */
        C_STATIONARY_SOLUTION = 2,                      /**< The solution corresponds to a Clarke stationary point. */
        M_STATIONARY_SOLUTION = 3,                      /**< The solution corresponds to a Mordukhovich stationary point. */
        S_STATIONARY_SOLUTION = 4                       /**< The solution corresponds to a strongly stationary point. */
    };

    enum printLevel {
        NONE = 0,                                       /**< No Output. */
        OUTER_LOOP_ITERATES = 1,                        /**< Print stats for each outer loop iterate. */
        INNER_LOOP_ITERATES = 2,                        /**< Print stats for each inner loop iterate. */
        SUBPROBLEM_SOLVER_ITERATES = 3                  /**< Print stats for each inner loop (and possibly output of subproblem solver). */
    };


    class Options {

        public:

            /** Default constructor. */
            Options( );


            /** Copy constructor (deep copy).
             *
             * @param rhs The object to be copied.
            */
            Options( const Options& rhs );


            /** Destructor. */
            ~Options( );


            /** Assignment operator.
             *
             * @param rhs The obejct from which to assign.
            */
            Options& operator=( const Options& rhs );


            /** Sets all options to default values. */
            void setToDefault( );


            /** Ensures the consistency of given options. */
            returnValue ensureConsistency( );

            double stationarityTolerance;               /**< Tolerance for 1-Norm of stationarity violation. */
            double complementarityTolerance;		    /**< Complementarity tolerance. */
            double initialComplementarityPenalty;	    /**< Start value for complementarity penalty term. */
            double complementarityPenaltyUpdate;	    /**< Factor for updating penaltised complementarity term. */

            bool solveZeroPenaltyFirst;                 /**< Flag indicating whether first QP should ignore penalization. */

            int maxOuterIterations;                     /**< Maximum number of outer iterations to be performed. */
            int maxInnerIterations;                     /**< Maximum number of inner iterations to be performed. */

            int relaxOptionsTolerance;                  /**< Number of active set changes until making subsolver options more percise. */

            printLevel printLvl;                        /**< Print level. */

        protected:
            void copy( const Options& rhs );        /**< Copy each property. */
    };

    class Utilities {
        public:
            // C = A*B.
            static void MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p);

            // C = A'*B
            static void TransponsedMatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p);

            // C = A'*B + B'*A
            static void MatrixSymmetrizationProduct(const double* const A, const double* const B, double* C, int m, int n);

            // d = A*b + c
            static void AffineLinearTransformation(const double alpha, const double* const A, const double* const b, const double* const c, double* d, int m, int n);

            // C = alpha*A + beta*B
            static void WeightedMatrixAdd(const double alpha, const double* const A, const double beta, const double* const B, double* C, int m, int n);

            // c = alpha*a + beta*b
            static void WeightedVectorAdd(const double alpha, const double* const a, const double beta, const double* const b, double* c, int m);

            // returns p' * Q * p
            static double QuadraticFormProduct(const double* const Q, const double* const p, int m);

            // returns a'*b
            static double DotProduct(const double* const a, const double* const b, int m);

            // returns 1-norm
            static double MaxAbs(const double* const a, int m);

            // Read integral data from file
            static returnValue readFromFile(int* data, int n, const char* datafilename);

            // Read float data from file
            static returnValue readFromFile(double* data, int n, const char* datafilename );

            // Read float data from file
            static returnValue writeToFile(double* data, int n, const char* datafilename );

            // Print a double valued matrix
            static void printMatrix(const double* const A, int m, int n, const char* const name);

            // Print an integer valued matrix
            static void printMatrix(const int* const A, int m, int n, const char* const name);

            // Printing bounds
            static void printStep(double* xk, double* pk, double* xk_new, double alpha, int nV);

            // Printing bounds
            static void printBounds(double* lb, double* xk, double* ub, int m);

            /** Transform a csc matrix to dense.
             *
             * @param sparse A sparse matrix.
             * @param full A target pointer for the full matrix (expected to have size m*n).
             * @param m Number of rows of `H_sparse` (in its dense representation).
             * @param n Number of columns of `H_sparse` (in its dense representation).
             *
             * @returns returnValue::SUCCESSFUL_RETURN, or returnValue::INDEX_OUT_OF_BOUNDSA if an index leads to invalid memory access of the dense array.
             */
            static returnValue csc_to_dns(const csc* const sparse, double* full, int m, int n);


            /** Transform a dense matrix to csc.
             *
             * @param full A dense double array.
             * @param m Number of rows of `full`.
             * @param n Number of columns of `full`.
             *
             * @returns A csc pointer to the sparse matrix.
             */
            static csc* dns_to_csc(const double* const full, int m, int n);


            // Methods below this line where taken from qpOASES implementation
            /** Returns the absolute value of a real number.
             * \return	Absolute value of a real number */
            inline double getAbs(double x);

            /** Checks if the absolute difference between x and y is less than TOL. */
            inline bool isEqual(double x, double y, double TOL);

            /** Checks if the absolute difference between x and y is less than the constant Utilities::ZERO. */
            inline bool isEqual(double x, double y);

            /** Checks if the absolute value of x is less than TOL. */
            inline bool isZero(double x, double TOL);

            /** Checks if the absolute value of x is less than the constant Utilities::ZERO. */
            inline bool isZero(double x);

            /** Returns the sign of a double x. */
            inline double getSign(double x);

            /** Returns the maximum of two integer values x and y. */
            inline int getMax(int x, int y);

            /** Returns the minimum of two integer values x and y. */
            inline int getMin(int x, int y);

            /** Returns the maximum of two real values x and y. */
            inline double getMax(double x, double y);

            /** Returns the minimum of two real values x and y. */
            inline double getMin(double x, double y);

            /** Returns the absolute value of x. */
            inline double getAbs(double x);

            /** Numerical value of machine precision (min eps, s.t. 1+eps > 1).
             *	Note: this value has to be positive! */
            #ifdef __USE_SINGLE_PRECISION__
            const double EPS = 1.193e-07f;
            #else
            const static double EPS = 2.221e-16;
            #endif /* __USE_SINGLE_PRECISION__ */

            /** Numerical value of zero (for situations in which it would be
             *	unreasonable to compare with 0.0).
            *	Note: this value has to be positive! */
            const static double ZERO = 1.0e-25;

            /** Numerical value of infinity (e.g. for non-existing bounds).
                Note: this value has to be positive! */
            const static double INFTY = 1.0e20;

            /** Maximum number of characters within a string.
             *	Note: this value should be at least 41! */
            const static uint MAX_STRING_LENGTH = 160;
    };

    class MessageHandler {
        public:
            static returnValue PrintMessage( returnValue ret );

            static algorithmStatus PrintSolution( algorithmStatus algoStat );

            static void PrintSolutionLine( );
    };
}

#include "Utilities.ipp"

#endif  // LCQPOASES_UTILITIES_HPP