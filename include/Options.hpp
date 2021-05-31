/*
 *	This file is part of LCQPanther.
 *
 *	LCQPanther -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	LCQPanther is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPanther is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPanther; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef LCQPanther_OPTIONS_HPP
#define LCQPanther_OPTIONS_HPP

#include "Utilities.hpp"

namespace LCQPanther {

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

            /** Get stationarity tolerance. */
            double getStationarityTolerance( );

            /** Set stationarity tolerance. */
            ReturnValue setStationarityTolerance( double val );

            /** Get complementarity tolerance. */
            double getComplementarityTolerance( );

            /** Set complementarity tolerance. */
            ReturnValue setComplementarityTolerance( double val );

            /** Get initial penalty parameter. */
            double getInitialPenaltyParameter( );

            /** Set complementarity tolerance. */
            ReturnValue setInitialPenaltyParameter( double val );

            /** Get penalty parameter update factor. */
            double getPenaltyUpdateFactor( );

            /** Set penalty parameter update factor. */
            ReturnValue setPenaltyUpdateFactor( double val );

            /** Get whether to solve for (complement.) unconstrained global minumum first. */
            bool getSolveZeroPenaltyFirst( );

            /** Set whether to solve for (complement.) unconstrained global minumum first. */
            ReturnValue setSolveZeroPenaltyFirst( bool val );

            /** Get maximum number of iterations. */
            int getMaxIterations( );

            /** Set maximum number of iterations. */
            ReturnValue setMaxIterations( int val );

            /** Get print level. */
            PrintLevel getPrintLevel( );

            /** Set print level. */
            ReturnValue setPrintLevel( PrintLevel val );

            /** Set print level (using an integer). */
            ReturnValue setPrintLevel( int val );

            /** Get QP solver. */
            QPSolver getQPSolver( );

            /** Set print level. */
            ReturnValue setQPSolver( QPSolver val );

            /** Set print level (using an integer). */
            ReturnValue setQPSolver( int val );

        protected:
            void copy( const Options& rhs );            /**< Copy each property. */

            double stationarityTolerance;               /**< Tolerance for 1-Norm of stationarity violation. */
            double complementarityTolerance;		    /**< Complementarity tolerance. */
            double initialPenaltyParameter;	            /**< Start value for complementarity penalty term. */
            double penaltyUpdateFactor;	                /**< Factor for updating penaltised complementarity term. */

            bool solveZeroPenaltyFirst;                 /**< Flag indicating whether first QP should ignore penalization. */

            int maxIterations;                          /**< Maximum number of iterations to be performed. */

            PrintLevel printLevel;                      /**< Print level. */

            QPSolver qpSolver;                          /**< The QP solver to be used. */
    };
}

#endif  // LCQPanther_UTILITIES_HPP
