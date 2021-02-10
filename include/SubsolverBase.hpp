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

#ifndef LCQPOASES_SUBSOLVERBASE_HPP
#define LCQPOASES_SUBSOLVERBASE_HPP

#include <Utilities.hpp>

namespace lcqpOASES {
    class SubsolverBase {
        public:
			/** Write solution to x. */
            virtual void getPrimalSolution( double* x ) = 0;

            /** Write solution to y. */
            virtual void getDualSolution( double* y ) = 0;

            /** Abstract method for solving the QP. */
            virtual returnValue solve(  bool initialSolve, int& iterations, 
                                        double* _g, 
                                        double* _lb, double* _ub,
                                        double* _lbA, double* _ubA  ) = 0;

    };
}

#endif  // LCQPOASES_SUBSOLVERBASE_HPP