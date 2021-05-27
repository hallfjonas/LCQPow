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

#include "Subsolver.hpp"
#include <cstring>

namespace LCQPanther {


    Subsolver::Subsolver( ) { }


    Subsolver::Subsolver(   int nV, int nC,
                            double* H, double* A )
    {
        qpSolver = QPSolver::QPOASES_DENSE;

        SubsolverQPOASES tmp(nV, nC, H, A);
        solverQPOASES = tmp;
    }


    Subsolver::Subsolver(   int nV, int nC,
                            csc* H, csc* A )
    {
        qpSolver = QPSolver::QPOASES_SPARSE;

        SubsolverQPOASES tmp(nV, nC, H, A );
        solverQPOASES = tmp;
    }


    Subsolver::Subsolver(   int nV, int nC,
                            csc* H, csc* A,
                            const double* g,
                            const double* l,
                            const double* u )
    {
        qpSolver = QPSolver::OSQP_SPARSE;

        SubsolverOSQP tmp(H, A, g, l, u);
        solverOSQP = tmp;
    }


    Subsolver::Subsolver(const Subsolver& rhs)
    {
        copy( rhs );
    }


    Subsolver& Subsolver::operator=(const Subsolver& rhs)
    {
        if ( this != &rhs )
            {
                copy( rhs );
            }

        return *this;
    }


    void Subsolver::getSolution( double* x, double* y )
    {
        if (qpSolver == QPSolver::QPOASES_DENSE || qpSolver == QPSolver::QPOASES_SPARSE) {
            solverQPOASES.getSolution( x, y );
        } else if (qpSolver == QPSolver::OSQP_SPARSE) {
            solverOSQP.getSolution( x, y );
        }
    }


    void Subsolver::setPrintLevel( PrintLevel printLevel )
    {
        if (qpSolver == QPSolver::QPOASES_DENSE || qpSolver == QPSolver::QPOASES_SPARSE) {
            if (printLevel < PrintLevel::SUBPROBLEM_SOLVER_ITERATES)
                optionsQPOASES.printLevel =  qpOASES::PrintLevel::PL_NONE;
            else
                optionsQPOASES.printLevel =  qpOASES::PrintLevel::PL_MEDIUM;

            solverQPOASES.setOptions( optionsQPOASES );
        } else if (qpSolver == QPSolver::OSQP_SPARSE) {
            solverOSQP.setPrintlevl(printLevel >= PrintLevel::SUBPROBLEM_SOLVER_ITERATES);
        }
    }


    ReturnValue Subsolver::solve(   bool initialSolve, int& iterations,
                                    const double* g,
                                    const double* lbA, const double* ubA,
                                    const double* x0, const double* y0,
                                    const double* lb, const double* ub)
    {
        ReturnValue ret = ReturnValue::SUCCESSFUL_RETURN;
        if (qpSolver == QPSolver::QPOASES_DENSE || qpSolver == QPSolver::QPOASES_SPARSE) {
            ret = solverQPOASES.solve( initialSolve, iterations, g, lbA, ubA, x0, y0, lb, ub );
        } else if (qpSolver == QPSolver::OSQP_SPARSE) {
            ret = solverOSQP.solve( initialSolve, iterations, g, lbA, ubA, x0, y0, lb, ub );
        } else {
            ret = INVALID_QPSOLVER;
        }

        return ret;
    }


    void Subsolver::copy(const Subsolver& rhs)
    {
        qpSolver = rhs.qpSolver;

        if (qpSolver == QPSolver::QPOASES_DENSE || qpSolver == QPSolver::QPOASES_SPARSE) {
            SubsolverQPOASES tmp( rhs.solverQPOASES );
            solverQPOASES = tmp;
        } else if (qpSolver == QPSolver::OSQP_SPARSE) {
            SubsolverOSQP tmp( rhs.solverOSQP );
            solverOSQP = tmp;
        }
    }
}