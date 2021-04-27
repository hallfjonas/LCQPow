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


    SubsolverOSQP::SubsolverOSQP( ) {
        data = 0;
        settings = 0;
        work = 0;
     }


    SubsolverOSQP::SubsolverOSQP(   const csc* const _H, const csc* const _A,
                                    const double* const _g,
                                    const double* const _l,
                                    const double* const _u)
    {
        // Store dimensions
        nVars = _H->n;
        nDuals = _A->m;

        // Allocate memory
        settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        data = (OSQPData *)c_malloc(sizeof(OSQPData));

        // Copy matrices
        H = copy_csc_mat(_H);
        A = copy_csc_mat(_A);
        g = new c_float[nVars];
        l = new c_float[nDuals];
        u = new c_float[nDuals];
        memcpy(g, _g, nVars*sizeof(c_float));
        memcpy(l, _l, nDuals*sizeof(c_float));
        memcpy(u, _u, nDuals*sizeof(c_float));

        // Fill data
        data->n = nVars;
        data->m = nDuals;
        data->P = H;
        data->A = A;
        data->q = g;
        data->l = l;
        data->u = u;

        // Define solver settings
        osqp_set_default_settings(settings);
        settings->eps_prim_inf = Utilities::ZERO;
        settings->verbose = false;

        // Setup workspace
        osqp_setup(&work, data, settings);
    }


    SubsolverOSQP::SubsolverOSQP(const SubsolverOSQP& rhs)
    {
        copy( rhs );
    }


    SubsolverOSQP::~SubsolverOSQP()
    {
        clear();
    }


    void SubsolverOSQP::clear()
    {
        osqp_cleanup(work);

        if (data != 0) {
            c_free(data);
            data = 0;
        }

        if (settings != 0) {
            c_free(settings);
            settings = 0;
        }

        if (H != 0) {
            csc_spfree(H);
            H = 0;
        }

        if (A != 0) {
            csc_spfree(A);
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


    SubsolverOSQP& SubsolverOSQP::operator=(const SubsolverOSQP& rhs)
    {
        if (this != &rhs) {
            copy( rhs );
        }

        return *this;
    }


    void SubsolverOSQP::setOptions( OSQPSettings* _settings )
    {
        settings = copy_settings(_settings);
    }


    void SubsolverOSQP::setPrintlevl( bool verbose )
    {
        if (settings != 0)
            settings->verbose = verbose;
    }


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

        iterations = work->info->iter;

        // Either pass error
        if (exitflag != 0)
            return returnValue::SUBPROBLEM_SOLVER_ERROR;

        // Or pass successful return
        return returnValue::SUCCESSFUL_RETURN;
    }


    void SubsolverOSQP::getSolution( double* x, double* y )
    {
        OSQPSolution *sol(work->solution);

        if (sol->x != 0) {
            memcpy(x, sol->x, nVars*(sizeof(double)));
        }

        // Copy duals with negative sign
        for (int i = 0; i < nDuals; i++) {
            y[i] = -sol->y[i];
        }
    }


    void SubsolverOSQP::copy(const SubsolverOSQP& rhs)
    {
        clear();

        nVars = rhs.nVars;
        nDuals = rhs.nDuals;

        H = copy_csc_mat(rhs.H);
        A = copy_csc_mat(rhs.A);
        g = new c_float[nVars];
        l = new c_float[nDuals];
        u = new c_float[nDuals];
        memcpy(g, rhs.g, nVars*sizeof(c_float));
        memcpy(l, rhs.l, nDuals*sizeof(c_float));
        memcpy(u, rhs.u, nDuals*sizeof(c_float));

        // Copy data
        data = (OSQPData *)c_malloc(sizeof(OSQPData));
        data->n = rhs.data->n;
        data->m = rhs.data->m;
        data->P = H;
        data->A = A;
        data->q = g;
        data->l = l;
        data->u = u;

        // Copy settings
        settings = copy_settings(rhs.settings);

        // Setup workspace
        osqp_setup(&work, data, settings);
    }
}