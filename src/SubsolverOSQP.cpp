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

#include <SubsolverOSQP.hpp>
#include "osqp.h"

namespace lcqpOASES {

    /*
     *   S u b s o l v e r O S Q P
     */
    SubsolverOSQP::SubsolverOSQP( ) { }


    /*
     *   S u b s o l v e r O  S Q P
     */
    SubsolverOSQP::SubsolverOSQP(   int nV, int nC,
                                    csc* _H, csc* _A )
    {
        settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        data = (OSQPData *)c_malloc(sizeof(OSQPData));

        H = _H;
        A = _A;

        // Populate data
        if (data) {
            data->n = nV;
            data->m = nC;
            data->P = H;
            data->A = A;
        }

        // Define solver settings as default
        if (settings) {
            osqp_set_default_settings(settings);
        }

        // Setup workspace
        osqp_setup(&work, data, settings);
    }

    /*
     *   S u b s o l v e r O  S Q P
     */
    SubsolverOSQP::SubsolverOSQP(const SubsolverOSQP& rhs)
    {
        copy( rhs );
    }


    /*
     *   ~ S u b s o l v e r O S Q P
     */
    SubsolverOSQP::~SubsolverOSQP()
    {
        // TODO: Destruct csc matrices
    }


    /*
     *   o p e r a t o r =
     */
    SubsolverOSQP& SubsolverOSQP::operator=(const SubsolverOSQP& rhs)
    {
        if (this != &rhs) {
            copy( rhs );
        }

        return *this;
    }


    /*
     *   s e t O p t i o n s
     */
    void SubsolverOSQP::setOptions( OSQPSettings* _settings )
    {
        settings = _settings;
    }


    /*
     *   s o l v e
     */
    returnValue SubsolverOSQP::solve(   bool initialSolve, int& iterations,
                                        const double* const _g,
                                        const double* const _lb, const double* const _ub,
                                        const double* const _lbA, const double* const _ubA,
                                        const double* const x0, const double* const y0 )
    {
        // Update linear cost and bounds
        osqp_update_lin_cost(work, _g);
        osqp_update_bounds(work, _lbA, _ubA);

        // Solve Problem
        int exitflag = osqp_solve(work);

        // Either pass error
        if (exitflag != 0)
            return returnValue::SUBPROBLEM_SOLVER_ERROR;

        // Or pass successful return
        return returnValue::SUCCESSFUL_RETURN;
    }


    /*
     *   g e t P r i m a l S o l u t i o n
     */
    void SubsolverOSQP::getPrimalSolution( double* x )
    {
        throw( MessageHandler::PrintMessage( returnValue::NOT_YET_IMPLEMENTED ));
    }


    /*
     *   g e t D u a l S o l u t i o n
     */
    void SubsolverOSQP::getDualSolution( double* y )
    {
        throw( MessageHandler::PrintMessage( returnValue::NOT_YET_IMPLEMENTED ));
    }


    /*
     *   c o p y
     */
    void SubsolverOSQP::copy(const SubsolverOSQP& rhs)
    {
        // Todo: Verify copy
        data = rhs.data;
        settings = rhs.settings;
        work = rhs.work;
    }
}