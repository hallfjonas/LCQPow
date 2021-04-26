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
    SubsolverOSQP::SubsolverOSQP( ) {
        data = 0;
        settings = 0;
        work = 0;
     }


    /*
     *   S u b s o l v e r O  S Q P
     */
    SubsolverOSQP::SubsolverOSQP(   int nV, int nC,
                                    csc* _H, csc* _A,
                                    const double* _g,
                                    const double* _l,
                                    const double* _u)
    {
        settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        data = (OSQPData *)c_malloc(sizeof(OSQPData));

        nVars = nV;
        nDuals = nC;

        H = csc_matrix(_H->m, _H->n, _H->nzmax, _H->x, _H->i, _H->p);
        A = csc_matrix(_A->m, _A->n, _A->nzmax, _A->x, _A->i, _A->p);

        // Conversion to OSQP data type c_float
        g = new c_float[nV];
        l = new c_float[nC];
        u = new c_float[nC];

        for (int i = 0; i < nV; i++) {
            g[i] = (c_float) _g[i];
        }

        for (int i = 0; i < nC; i++) {
            l[i] = (c_float) _l[i];
            u[i] = (c_float) _u[i];
        }

		Utilities::printMatrix(H, "Full H");
		Utilities::printMatrix(A, "Full A");

        // Populate data
        if (data) {
            data->n = nV;
            data->m = nC;
            data->P = H;
            data->A = A;
            data->q = g;
            data->l = l;
            data->u = u;
        }

        // Define solver settings as default
        if (settings) {
            osqp_set_default_settings(settings);
            settings->eps_prim_inf = Utilities::ZERO;
            settings->verbose = false;
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
        if (settings != 0) {
            c_free(settings);
            settings = 0;
        }

        if (data != 0) {
            c_free(data);
            data = 0;
        }

        // TODO: FIGURE OUT WHEN AND HOW TO CALL OSQP_CLEANUP
        printf("NOT CLEANING UP osqp work!!!\n");

        if (H != 0) {
            c_free(H);
            H = 0;
        }

        if (A != 0) {
            c_free(A);
            A = 0;
        }

        if (g != 0) {
            delete[] g;
            g = 0;
        }

        if (l != 0) {
            delete[] l;
            l = 0;
        }

        if (u != 0) {
            delete[] u;
            u = 0;
        }
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
     *   s e t O p t i o n s
     */
    void SubsolverOSQP::setPrintlevl( bool verbose )
    {
        if (settings != 0)
            settings->verbose = verbose;
    }


    /*
     *   s o l v e
     */
    returnValue SubsolverOSQP::solve(   bool initialSolve, int& iterations,
                                        const double* const _g,
                                        const double* const _lbA, const double* const _ubA,
                                        const double* const x0, const double* const y0,
                                        const double* const _lb, const double* const _ub )
    {
        // Make sure that lb and ub are null pointers, as OSQP does not handle box constraints
        if (_lb != 0 || _ub != 0) {
            return MessageHandler::PrintMessage( returnValue::INVALID_OSQP_BOX_CONSTRAINTS );
        }

        // Update linear cost and bounds
        osqp_update_lin_cost(work, _g);
        osqp_update_bounds(work, _lbA, _ubA);

        int exitflag;

        // Solve Problem
        exitflag = osqp_solve(work);

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
        OSQPSolution *sol(work->solution);

        if (sol->x != 0) {
            memcpy(x, sol->x, nVars*(sizeof(double)));
        }
    }


    /*
     *   g e t D u a l S o l u t i o n
     */
    void SubsolverOSQP::getDualSolution( double* y )
    {
        OSQPSolution *sol(work->solution);

        if (sol->y != 0) {
            memcpy(y, sol->y, nDuals*(sizeof(double)));
        }
    }


    /*
     *   c o p y
     */
    void SubsolverOSQP::copy(const SubsolverOSQP& rhs)
    {
        nVars = rhs.nVars;
        nDuals = rhs.nDuals;

        if (rhs.H != 0)
            H = csc_matrix(rhs.H->m, rhs.H->n, rhs.H->nzmax, rhs.H->x, rhs.H->i, rhs.H->p);

        if (rhs.A != 0)
            A = csc_matrix(rhs.A->m, rhs.A->n, rhs.A->nzmax, rhs.A->x, rhs.A->i, rhs.A->p);

        if (rhs.data) {
            data = (OSQPData *)c_malloc(sizeof(OSQPData));
            memcpy(data, rhs.data, sizeof(OSQPData));
        }

        if (rhs.settings) {
            settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
            memcpy(settings, rhs.settings, sizeof(OSQPSettings));
        }

        if (rhs.work) {
            work = (OSQPWorkspace *)c_malloc(sizeof(OSQPWorkspace));
            memcpy(work, rhs.work, sizeof(OSQPWorkspace));
        }
    }
}