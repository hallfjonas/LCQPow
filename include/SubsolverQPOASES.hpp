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

#ifndef LCQPanther_SUBSOLVERQPOASES_HPP
#define LCQPanther_SUBSOLVERQPOASES_HPP

#include "SubsolverBase.hpp"
#include <qpOASES.hpp>

namespace LCQPanther {
    class SubsolverQPOASES : public SubsolverBase {
        public:
			/** Default constructor. */
			SubsolverQPOASES( );

            /** Constructor. */
            SubsolverQPOASES(   int nV,
                                int nC,
                                double* H,
                                double* A);

            /** Constructor. */
            SubsolverQPOASES(   int nV,
                                int nC,
                                csc* H,
                                csc* A);

            /** Copy constructor. */
            SubsolverQPOASES(const SubsolverQPOASES& rhs);

            /** Destructor. */
            ~SubsolverQPOASES();

            /** Assignment operator (deep copy). */
            virtual SubsolverQPOASES& operator=(const SubsolverQPOASES& rhs);

            /** Setting the user options. */
            void setOptions( qpOASES::Options options );

            /** Run qpOASES solver. */
            ReturnValue solve(  bool initialSolve, int& iterations,
                                const double* const _g,
                                const double* const _lbA,
                                const double* const _ubA,
                                const double* const x0 = 0,
                                const double* const y0 = 0,
                                const double* const _lb = 0,
                                const double* const _ub = 0 );

            /** Write solution to x. */
            void getSolution( double* x, double* y );

        protected:
            /** Copies all members from given rhs object. */
            void copy(const SubsolverQPOASES& rhs);

        private:
            int nV;
            int nC;

            // Dense types
            double* H = NULL;
            double* A = NULL;

            bool isSparse = false;

            // Sparse types (must be of generic types for duplicate method)
            qpOASES::SymSparseMat* H_sparse = NULL;
            qpOASES::SparseMatrix* A_sparse = NULL;

            double* H_x = NULL;
            int* H_i = NULL;
            int* H_p = NULL;

            double* A_x = NULL;
            int* A_i = NULL;
            int* A_p = NULL;

            qpOASES::QProblem qp;
    };
}

#endif  // LCQPanther_SUBSOLVERQPOASES_HPP