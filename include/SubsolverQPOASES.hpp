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

#ifndef LCQPOASES_SUBSOLVERQPOASES_HPP
#define LCQPOASES_SUBSOLVERQPOASES_HPP

#include "SubsolverBase.hpp"
#include <qpOASES.hpp>

namespace lcqpOASES {
    class SubsolverQPOASES : public SubsolverBase {
        public:
			/** Default constructor. */
			SubsolverQPOASES( );

            /** Constructor. */
            SubsolverQPOASES(   int nV,
                                int nC,
                                double* H,
                                double* A);

            /** Copy constructor. */
            SubsolverQPOASES(const SubsolverQPOASES& rhs);

            /** Destructor. */
            ~SubsolverQPOASES();

            /** Assignment operator (deep copy). */
            virtual SubsolverQPOASES& operator=(const SubsolverQPOASES& rhs);

            /** Setting the user options. */
            void setOptions( qpOASES::Options options );

            /** Run qpOASES solver. */
            returnValue solve(  bool initialSolve, int& iterations,
                                const double* const _g,
                                const double* const _lbA,
                                const double* const _ubA,
                                const double* const x0 = 0,
                                const double* const y0 = 0,
                                const double* const _lb = 0,
                                const double* const _ub = 0 );

            /** Write solution to x. */
            void getPrimalSolution( double* x );

            /** Write solution to y. */
            void getDualSolution( double* y );

        protected:
            /** Copies all members from given rhs object. */
            void copy(const SubsolverQPOASES& rhs);

        private:
            int nV;
            int nC;

            double* H;
            double* A;

            qpOASES::QProblem qp;
    };
}

#endif  // LCQPOASES_SUBSOLVERQPOASES_HPP