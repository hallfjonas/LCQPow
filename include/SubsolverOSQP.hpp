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

#ifndef LCQPOASES_SUBSOLVEROSQP_HPP
#define LCQPOASES_SUBSOLVEROSQP_HPP

#include "SubsolverBase.hpp"
#include <osqp.h>

namespace lcqpOASES {
    class SubsolverOSQP : public SubsolverBase {
        public:
			/** Default constructor. */
			SubsolverOSQP( );

            SubsolverOSQP(  int _nV,
                            int _nC,
                            csc* _H,
                            csc* _A,
                            const double* g,
                            const double* lbA,
                            const double* ubA
                            );

            /** Copy constructor. */
            SubsolverOSQP(const SubsolverOSQP& rhs);

            /** Destructor. */
            ~SubsolverOSQP( );

            /** Assignment operator (deep copy). */
            virtual SubsolverOSQP& operator=(const SubsolverOSQP& rhs);

            /** Set OSQP settings. */
            void setOptions( OSQPSettings* settings );

            /** Run OSQP solver. */
            returnValue solve(  bool initialSolve, int& iterations,
                                const double* const _g,
                                const double* const _lb, const double* const _ub,
                                const double* const _lbA, const double* const _ubA,
                                const double* const x0 = 0, const double* const y0 = 0);

            /** Write solution to x. */
            void getPrimalSolution( double* x );

            /** Write solution to y. */
            void getDualSolution( double* y );

        protected:
            /** Copies all members from given rhs object. */
            void copy(const SubsolverOSQP& rhs);

        private:
            csc* H;
            csc* A;

            OSQPWorkspace *work;
            OSQPSettings *settings;
            OSQPData *data;
    };
}

#endif  // LCQPOASES_SUBSOLVEROSQP_HPP