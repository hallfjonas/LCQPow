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

#ifndef LCQPOASES_SUBSOLVER_HPP
#define LCQPOASES_SUBSOLVER_HPP

#include "SubsolverQPOASES.hpp"
#include "SubsolverOSQP.hpp"

namespace lcqpOASES {

    enum QPSubproblemSolver {
        QPOASES = 0,                                /**< qpOASES. */
        OSQP = 1
    };

    class Subsolver {
        public:
			/** Default constructor. */
			Subsolver( );

            /** Constructor for dense matrices. */
            Subsolver(  int nV,
                        int nC,
                        double* H,
                        double* A );

            /** Constructor for sparse matrices. */
            Subsolver(  int nV,
                        int nC,
                        csc* H,
                        csc* A,
                        const double* g,
                        const double* l,
                        const double* u);

            /** Copy constructor. */
            Subsolver(const Subsolver& rhs);

            /** Assignment operator (deep copy). */
            virtual Subsolver& operator=(const Subsolver& rhs);

            /** Write solution to x and y. */
            void getSolution( double* x, double* y );

            /** Options of subproblem solver. */
            void setPrintLevel( printLevel printlvl );

            /** Abstract method for solving the QP. */
            returnValue solve(  bool initialSolve, int& iterations,
                                const double* g,
                                const double* lbA, const double* ubA,
                                const double* x0 = 0, const double* y0 = 0,
                                const double* lb = 0, const double* ub = 0);

        protected:
            /** Copies all members from given rhs object. */
            void copy(const Subsolver& rhs);


        private:
            // The solver type
            QPSubproblemSolver qpSolver;        	/**< Inidicating which qpSolver to use. */

            // The different solvers
        	SubsolverQPOASES solverQPOASES;			/**< When using qpOASES. */
			SubsolverOSQP solverOSQP;				/**< When using OSQP. */

            // Options and settings
            qpOASES::Options optionsQPOASES;        /**< Options for the qpOASES solver. */
    };
}

#endif  // LCQPOASES_SUBSOLVER_HPP