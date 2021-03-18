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

#include <Subsolver.hpp>
#include <cstring>

namespace lcqpOASES {

    Subsolver::Subsolver( ) { }

    /*
     *   S u b s o l v e r (QPOASES)
     */
    Subsolver::Subsolver(   int nV, int nC,
                            double* H, double* A )
    {
        qpSolver = QPSubproblemSolver::QPOASES;
        SubsolverQPOASES tmp(nV, nC, H, A);
        solverQPOASES = tmp;
    }

    /*
     *   S u b s o l v e r (OSQP)
     */
    Subsolver::Subsolver(   int nV, int nC,
                            csc* H, csc* A )
    {
        qpSolver = QPSubproblemSolver::OSQP;
        solverOSQP = SubsolverOSQP(nV, nC, H, A);
    }

    /** Copy constructor. */
    Subsolver::Subsolver(const Subsolver& rhs)
    {
        copy( rhs );
    }

    /** Assignment operator (deep copy). */
    Subsolver& Subsolver::operator=(const Subsolver& rhs)
    {
        if ( this != &rhs )
            {
                copy( rhs );
            }

        return *this;
    }


    /** Write solution to x. */
    void Subsolver::getPrimalSolution( double* x )
    {
        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
                solverQPOASES.getPrimalSolution( x );
                return;

            case QPSubproblemSolver::OSQP:
                solverQPOASES.getPrimalSolution( x );
                return;
        }
    }

    /** Write solution to y. */
    void Subsolver::getDualSolution( double* y )
    {
        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
                solverQPOASES.getDualSolution( y );
                return;

            case QPSubproblemSolver::OSQP:
                solverQPOASES.getDualSolution( y );
                return;
        }
    }


    /*
     *   s e t o p t i o n s
     */
    void Subsolver::setPrintLevel( printLevel printlvl )
    {
        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
            {
                if (printlvl < printLevel::SUBPROBLEM_SOLVER_ITERATES)
                    optionsQPOASES.printLevel =  qpOASES::PrintLevel::PL_NONE;
                else
                    optionsQPOASES.printLevel =  qpOASES::PrintLevel::PL_MEDIUM;

                solverQPOASES.setOptions( optionsQPOASES );
                break;
            }

            case QPSubproblemSolver::OSQP:
            {
                MessageHandler::PrintMessage( NOT_YET_IMPLEMENTED );
                return;
            }
        }
    }


    /*
     *   s w i t c h T o R e l a x e d O p t i o n s
     */
    void Subsolver::switchToRelaxedOptions( )
    {
        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
            {
                qpOASES::PrintLevel pl = optionsQPOASES.printLevel;
                optionsQPOASES.setToFast( );
                optionsQPOASES.printLevel = pl;
                solverQPOASES.setOptions( optionsQPOASES );
                break;
            }

            case QPSubproblemSolver::OSQP:
            {
                MessageHandler::PrintMessage( NOT_YET_IMPLEMENTED );
                break;
            }
        }
    }


    /*
     *   s w i t c h T o S t r i c t O p t i o n s
     */
    void Subsolver::switchToStrictOptions( )
    {
        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
            {
                qpOASES::PrintLevel pl = optionsQPOASES.printLevel;
                optionsQPOASES.setToDefault( );
                optionsQPOASES.printLevel = pl;
                solverQPOASES.setOptions( optionsQPOASES );
                break;
            }

            case QPSubproblemSolver::OSQP:
            {
                MessageHandler::PrintMessage( NOT_YET_IMPLEMENTED );
                break;
            }
        }
    }


    /*
     *   s o l v e
     */
    returnValue Subsolver::solve(   bool initialSolve, int& iterations,
                                    double* g,
                                    double* lb, double* ub,
                                    double* lbA, double* ubA,
                                    double* x0, double* y0 )
    {
        returnValue ret = returnValue::SUCCESSFUL_RETURN;
        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
                ret = solverQPOASES.solve( initialSolve, iterations, g, lb, ub, lbA, ubA, x0, y0 );
                break;

            case QPSubproblemSolver::OSQP:
                ret = solverOSQP.solve( initialSolve, iterations, g, lb, ub, lbA, ubA, x0, y0 );
                break;
        }

        return ret;
    }


    /*
     *   c o p y
     */
    void Subsolver::copy(const Subsolver& rhs)
    {
        qpSolver = rhs.qpSolver;

        switch (qpSolver) {
            case QPSubproblemSolver::QPOASES:
                solverQPOASES = SubsolverQPOASES( rhs.solverQPOASES );
                return;

            case QPSubproblemSolver::OSQP:
                solverOSQP = SubsolverOSQP( rhs.solverOSQP );
                return;
        }
    }

}