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


#ifndef LCQPow_OUTPUTSTATISTICS_HPP
#define LCQPow_OUTPUTSTATISTICS_HPP

#include "Utilities.hpp"

namespace LCQPow {

    class OutputStatistics {
        public:
            /** Default constructor. */
            OutputStatistics( );

            /** Assignment operator.
             *
             * @param rhs The obejct from which to assign.
            */
            OutputStatistics& operator=( const OutputStatistics& rhs );


            /** Resets the statistics. */
            void reset( );

            /** Update total iteration counter.
             *
             * @return Success or specifies the invalid argument.
            */
            ReturnValue updateIterTotal( int delta_iter );

            /** Update total outer iteration counter.
             *
             * @return Success or specifies the invalid argument.
            */
            ReturnValue updateIterOuter( int delta_iter );

            /** Update total number of working set changes counter.
             *
             * @return Success or specifies the invalid argument.
            */
            ReturnValue updateSubproblemIter( int delta_iter );

            /** Update rho at solution.
             *
             * @return Success or specifies the invalid argument.
            */
            ReturnValue updateRhoOpt( double _rho );

            /** Update the solution status.
             *
             * @return Success or specifies the invalid argument.
            */
            ReturnValue updateSolutionStatus( AlgorithmStatus _status );

            /** Get the total number of iterations. */
            int getIterTotal( ) const;

            /** Get the total number of outer iterations. */
            int getIterOuter( ) const;

            /** Get the total number of subproblem iterations. */
            int getSubproblemIter( ) const;

            /** Get the penalty parameter at the optimal solution (if found). */
            double getRhoOpt( ) const;

            /** Get the solution status (if solved it will return the stationarity type). */
            AlgorithmStatus getSolutionStatus( ) const;

        private:
            int iter_total = 0;                                 /**< Total number of iterations, i.e., total number of inner iterations. */
            int iter_outer = 0;                                 /**< Total number of outer iterations, i.e., number of penalty updates. */
            int subproblem_iter = 0;                            /**< Total number of subsolver iterations (qpOASES: active set changes, OSQP: ??). */
            double rho_opt = 0.0;                               /**< Value of penalty parameter at the final iterate. */
            AlgorithmStatus status = PROBLEM_NOT_SOLVED;        /**< Status of the solver. This is set to the solution type on success. */
    };
}

#endif  // LCQPow_OUTPUTSTATISTICS_HPP