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

#include <SubsolverQPOASES.hpp>
#include <qpOASES.hpp>

namespace LCQPanther {

    /*
     *   S u b s o l v e r O S Q P
     */
    SubsolverQPOASES::SubsolverQPOASES( ) { }


    /*
     *   S u b s o l v e r Q P O A S E S
     */
    SubsolverQPOASES::SubsolverQPOASES( int _nV, int _nC,
                                        double* _H, double* _A)
    {
        nV = _nV;
        nC = _nC;

        isSparse = false;
        qp = qpOASES::QProblem(nV, nC);

        H = new double[nV*nV];
        A = new double[nC*nV];

        memcpy(H, _H, (size_t)(nV*nV)*sizeof(double));
        memcpy(A, _A, (size_t)(nC*nV)*sizeof(double));
    }


    SubsolverQPOASES::SubsolverQPOASES( int _nV, int _nC,
                                        csc* _H, csc* _A)
    {
        nV = _nV;
        nC = _nC;

        isSparse = true;

        // Use this qpOASES flag to identify a sparse solver
        #ifndef SOLVER_NONE
            useSchur = true;
        #else
            useSchur = false;
        #endif

        if (useSchur) {
            qpSchur = qpOASES::SQProblemSchur(nV, nC);
        } else {
            qp = qpOASES::QProblem(nV, nC);
        }

        if (H_sparse != NULL) {
            free(H_sparse);
            H_sparse = NULL;
        }

        if (A_sparse != NULL) {
            free(A_sparse);
            A_sparse = NULL;
        }

        H_i = new int[_H->p[_nV]];
        H_x = new double[_H->p[_nV]];
        H_p = new int[_nV+1];

        A_i = new int[_A->p[_nV]];
        A_x = new double[_A->p[_nV]];
        A_p = new int[_nV+1];

        memcpy(H_p, _H->p, (size_t)(nV+1)*sizeof(int));
        memcpy(H_i, _H->i, (size_t)(_H->p[nV])*sizeof(int));
        memcpy(H_x, _H->x, (size_t)(_H->p[nV])*sizeof(double));
        memcpy(A_p, _A->p, (size_t)(nV+1)*sizeof(int));
        memcpy(A_i, _A->i, (size_t)(_A->p[nV])*sizeof(int));
        memcpy(A_x, _A->x, (size_t)(_A->p[nV])*sizeof(double));

        H_sparse = new qpOASES::SymSparseMat(nV, nV, H_i, H_p, H_x);
        A_sparse = new qpOASES::SparseMatrix(nC, nV, A_i, A_p, A_x);
        H_sparse->createDiagInfo();
        A_sparse->createDiagInfo();
    }

    /*
     *   S u b s o l v e r O  S Q P
     */
    SubsolverQPOASES::SubsolverQPOASES(const SubsolverQPOASES& rhs)
    {
        copy( rhs );
    }

    /*
     *   ~ S u b s o l v e r Q P O A S E S
     */
    SubsolverQPOASES::~SubsolverQPOASES()
    {
        if (H != 0) {
            delete[] H;
            H = NULL;
        }

        if (A != 0) {
            delete[] A;
            A = NULL;
        }

        if (H_i != 0) {
            delete[] H_i;
        }

        if (H_p != 0) {
            delete[] H_p;
        }

        if (H_x != 0) {
            delete[] H_x;
        }

        if (A_i != 0) {
            delete[] A_i;
        }

        if (A_p != 0) {
            delete[] A_p;
        }

        if (A_x != 0) {
            delete[] A_x;
        }

        if (H_sparse != 0) {
            delete H_sparse;
            H_sparse = 0;
        }

        if (A_sparse != 0) {
            delete A_sparse;
            A_sparse = 0;
        }
    }

    /*
     *   o p e r a t o r =
     */
    SubsolverQPOASES& SubsolverQPOASES::operator=(const SubsolverQPOASES& rhs)
    {
        if (this != &rhs) {
            copy( rhs );
        }

        return *this;
    }


    /*
     *   s e t O p t i o n s
     */
    void SubsolverQPOASES::setOptions( qpOASES::Options options )
    {
        if (useSchur) {
            qpSchur.setOptions( options );
        } else {
            qp.setOptions( options );
        }
    }


    /*
     *   s o l v e
     */
    ReturnValue SubsolverQPOASES::solve(    bool initialSolve, int& iterations,
                                            const double* const g,
                                            const double* const lbA, const double* const ubA,
                                            const double* const x0, const double* const y0,
                                            const double* const lb, const double* const ub )
    {
        qpOASES::returnValue ret;

        qpOASES::int_t nwsr = 1000000;

        if (initialSolve) {
            if (isSparse) {
                if (useSchur) {
                    ret = qpSchur.init(H_sparse, g, A_sparse, lb, ub, lbA, ubA, nwsr, (double*)0, x0, y0);
                } else {
                    ret = qp.init(H_sparse, g, A_sparse, lb, ub, lbA, ubA, nwsr, (double*)0, x0, y0);
                }
            } else {
                ret = qp.init(H, g, A, lb, ub, lbA, ubA, nwsr, (double*)0, x0, y0);
            }
        } else {
            if (useSchur) {
                ret = qpSchur.hotstart(g, lb, ub, lbA, ubA, nwsr);
            } else {
                ret = qp.hotstart(g, lb, ub, lbA, ubA, nwsr);
            }
        }

        iterations = (int)(nwsr);

        if (ret != qpOASES::returnValue::SUCCESSFUL_RETURN)
            return ReturnValue::SUBPROBLEM_SOLVER_ERROR;

        return ReturnValue::SUCCESSFUL_RETURN;
    }


    /*
     *   g e t P r i m a l S o l u t i o n
     */
    void SubsolverQPOASES::getSolution( double* x, double* y )
    {
        if (useSchur) {
            qpSchur.getPrimalSolution( x );
            qpSchur.getDualSolution( y );
        } else {
            qp.getPrimalSolution( x );
            qp.getDualSolution( y );
        }
    }


    /*
     *   c o p y
     */
    void SubsolverQPOASES::copy(const SubsolverQPOASES& rhs)
    {
        nV = rhs.nV;
        nC = rhs.nC;

        isSparse = rhs.isSparse;
        useSchur = rhs.useSchur;

        if (isSparse) {
            H_i = new int[rhs.H_p[nV]];
            H_x = new double[rhs.H_p[nV]];
            H_p = new int[nV+1];

            A_i = new int[rhs.A_p[nV]];
            A_x = new double[rhs.A_p[nV]];
            A_p = new int[nV+1];

            memcpy(H_p, rhs.H_p, (size_t)(nV+1)*sizeof(int));
            memcpy(H_i, rhs.H_i, (size_t)(rhs.H_p[nV])*sizeof(int));
            memcpy(H_x, rhs.H_x, (size_t)(rhs.H_p[nV])*sizeof(double));
            memcpy(A_p, rhs.A_p, (size_t)(nV+1)*sizeof(int));
            memcpy(A_i, rhs.A_i, (size_t)(rhs.A_p[nV])*sizeof(int));
            memcpy(A_x, rhs.A_x, (size_t)(rhs.A_p[nV])*sizeof(double));

            H_sparse = new qpOASES::SymSparseMat(nV, nV, H_i, H_p, H_x);
            A_sparse = new qpOASES::SparseMatrix(nC, nV, A_i, A_p, A_x);
            H_sparse->createDiagInfo();
            A_sparse->createDiagInfo();
        } else {
            H = new double[nV*nV];
            A = new double[nC*nV];

            memcpy(H, rhs.H, (size_t)(nV*nV)*sizeof(double));
            memcpy(A, rhs.A, (size_t)(nC*nV)*sizeof(double));
        }

        if (useSchur) {
            qpSchur = rhs.qpSchur;
        } else {
            qp = rhs.qp;
        }
    }
}