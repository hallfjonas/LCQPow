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

#ifndef LCQPanther_SUBSOLVERBASE_HPP
#define LCQPanther_SUBSOLVERBASE_HPP

#include "Utilities.hpp"

namespace LCQPanther {
    class SubsolverBase {
        public:
			/** Write solution to x. */
            virtual void getSolution( double* x, double* y ) = 0;

            /** Abstract method for solving the QP. */
            virtual ReturnValue solve(  bool initialSolve, int& iterations,
                                        const double* const _g,
                                        const double* const _lb, const double* const _ub,
                                        const double* const _lbA, const double* const _ubA,
                                        const double* const x0 = 0, const double* const y0 = 0) = 0;

    };
}

#endif  // LCQPanther_SUBSOLVERBASE_HPP