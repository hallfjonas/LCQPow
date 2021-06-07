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


#include "LCQProblem.hpp"
#include "Utilities.hpp"
#include "MessageHandler.hpp"
#include "SubsolverQPOASES.hpp"
#include "SubsolverOSQP.hpp"

#include <iostream>
#include <string>
#include <math.h>
#include <qpOASES.hpp>

using qpOASES::QProblem;

namespace LCQPanther {

	LCQProblem::LCQProblem( ) { }


	LCQProblem::LCQProblem( int _nV, int _nC, int _nComp )
	{
		/* consistency checks */
		if ( _nV <= 0 )
		{
			MessageHandler::PrintMessage( INVALID_NUMBER_OF_OPTIM_VARS );
			return;
		}

		if ( _nComp <= 0 )
		{
			MessageHandler::PrintMessage( INVALID_NUMBER_OF_OPTIM_VARS );
			return;
		}

		if ( _nC < 0 )
		{
			_nC = 0;
			MessageHandler::PrintMessage( INVALID_NUMBER_OF_CONSTRAINT_VARS );
			return;
		}

		nV = _nV;
		nC = _nC;
		nComp = _nComp;

		// Allocate auxiliar vectors
		Qk = new double[nV*nV]();
		gk = new double[nV]();
		xnew = new double[nV]();
		yk_A = new double[nC + 2*nComp]();
		pk = new double[nV]();
		statk = new double[nV]();
		constr_statk = new double[nV]();
		box_statk = new double[nV]();
		lk_tmp = new double[nV]();
	}


	LCQProblem::~LCQProblem( ) {
		clear();
	}


	ReturnValue LCQProblem::loadLCQP(	const double* const _H, const double* const _g,
										const double* const _S1, const double* const _S2,
										const double* const _A, const double* const _lbA, const double* const _ubA,
										const double* const _lb, const double* const _ub,
										const double* const _x0, const double* const _y0
										)
	{
		ReturnValue ret;

		if ( nV <= 0 || nComp <= 0 )
            return( MessageHandler::PrintMessage(ReturnValue::LCQPOBJECT_NOT_SETUP) );

		ret = setH( _H );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		// Only copy lb and ub to temporary variables
		// Build them once we know what solver is used
		if (_lb != 0) {
			lb_tmp = new double[nV];
			memcpy(lb_tmp, _lb, (size_t) nV*sizeof(double));
		}

		if (_ub != 0) {
			ub_tmp = new double[nV];
			memcpy(ub_tmp, _ub, (size_t) nV*sizeof(double));
		}

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setConstraints( _S1, _S2, _A, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		sparseSolver = false;

		return ReturnValue::SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::loadLCQP(	const char* const H_file, const char* const g_file,
										const char* const S1_file, const char* const S2_file,
										const char* const A_file, const char* const lbA_file, const char* const ubA_file,
										const char* const lb_file, const char* const ub_file,
										const char* const x0_file, const char* const y0_file
										)
	{
		ReturnValue ret;

		double* _H = new double[nV*nV];
		ret = Utilities::readFromFile( _H, nV*nV, H_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _H;
			return MessageHandler::PrintMessage( ret );
		}

		double* _g = new double[nV];
		ret = Utilities::readFromFile( _g, nV, g_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _g;
			return MessageHandler::PrintMessage( ret );
		}

		double* _S1 = new double[nComp*nV];
		ret = Utilities::readFromFile( _S1, nComp*nV, S1_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _S1;
			return MessageHandler::PrintMessage( ret );
		}

		double* _S2 = new double[nComp*nV];
		ret = Utilities::readFromFile( _S2, nComp*nV, S2_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _S2;
			return MessageHandler::PrintMessage( ret );
		}

		double* _A = NULL;
		if (A_file != 0) {
			_A = new double[nC*nV];
			ret = Utilities::readFromFile( _A, nC*nV, A_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _A;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _lbA = NULL;
		if (lbA_file != 0) {
			_lbA = new double[nC];
			ret = Utilities::readFromFile( _lbA, nC, lbA_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _lbA;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _ubA = NULL;
		if (ubA_file != 0) {
			_ubA = new double[nC];
			ret = Utilities::readFromFile( _ubA, nC, ubA_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _ubA;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _lb = NULL;
		if (lb_file != 0) {
			_lb = new double[nV];
			ret = Utilities::readFromFile( _lb, nV, lb_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _lb;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _ub = NULL;
		if (ub_file != 0) {
			_ub = new double[nV];
			ret = Utilities::readFromFile( _ub, nV, ub_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _ub;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _x0 = NULL;
		if (x0_file != 0) {
			_x0 = new double[nV];
			ret = Utilities::readFromFile( _x0, nV, x0_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _x0;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _y0 = NULL;
		if (y0_file != 0) {
			_y0 = new double[nC + 2*nComp];
			ret = Utilities::readFromFile( _y0, nC + 2*nComp, y0_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _y0;
				return MessageHandler::PrintMessage( ret );
			}
		}

		// Fill vaues
		ret = setH( _H );
		delete[] _H;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setG( _g );
		delete[] _g;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		// Only copy lb and ub to temporary variables
		// Build them once we know what solver is used
		if (_lb != 0) {
			lb_tmp = new double[nV];
			memcpy(lb_tmp, _lb, (size_t) nV*sizeof(double));
			delete[] _lb;
		}

		if (_ub != 0) {
			ub_tmp = new double[nV];
			memcpy(ub_tmp, _ub, (size_t) nV*sizeof(double));
			delete[] _ub;
		}

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setConstraints( _S1, _S2, _A, _lbA, _ubA );

		delete[] _S1; delete[] _S2;

		if (_A != 0)
			delete[] _A;

		if (_lbA != 0)
			delete[] _lbA;

		if (_ubA != 0)
			delete[] _ubA;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setInitialGuess( _x0, _y0 );

		if (_x0 != 0) delete[] _x0;
		if (_y0 != 0) delete[] _y0;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		sparseSolver = false;

		return ReturnValue::SUCCESSFUL_RETURN;
	}




	ReturnValue LCQProblem::loadLCQP(	double* _H_data, int* _H_i, int* _H_p, double* _g,
										double* _S1_data, int* _S1_i, int* _S1_p,
										double* _S2_data, int* _S2_i, int* _S2_p,
										double* _A_data, int* _A_i, int* _A_p,
										double* _lbA, double* _ubA,	double* _lb, double* _ub,
										double* _x0, double* _y0
										)
	{
		ReturnValue ret;

		ret = setH( _H_data, _H_i, _H_p );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setConstraints( _S1_data, _S1_i, _S1_p, _S2_data, _S2_i, _S2_p, _A_data, _A_i, _A_p, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		// Only copy lb and ub to temporary variables
		// Build them once we know what solver is used
		if (_lb != 0) {
			lb_tmp = new double[nV];
			memcpy(lb_tmp, _lb, (size_t) nV*sizeof(double));
		}

		if (_ub != 0) {
			ub_tmp = new double[nV];
			memcpy(ub_tmp, _ub, (size_t) nV*sizeof(double));
		}

		sparseSolver = true;

		return ReturnValue::SUCCESSFUL_RETURN;
	}

	ReturnValue LCQProblem::runSolver( )
	{
		// Initialize variables
		ReturnValue ret = initializeSolver();
		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage(ret);

		// Initialization strategy
		if (options.getSolveZeroPenaltyFirst()) {
			memcpy(gk, g, (size_t)nV*sizeof(double));

			ret = solveQPSubproblem( true );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage(ret);
			}

			updateLinearization();

			ret = solveQPSubproblem( false );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage(ret);
			}

		} else {
			updateLinearization();

			ret = solveQPSubproblem( true );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage(ret);
			}
		}

		// Outer and inner loop in one
		while ( true ) {

			// Update xk, gk, Qk, stationarity
			updateStep( );

			// Print iteration
			printIteration( );

			// Terminate, update pen, or continue inner loop
			if (stationarityCheck()) {
				if (complementarityCheck()) {
					// Switch from penalized to LCQP duals
					transformDuals();

					// Determine C-,M-,S-Stationarity
					determineStationarityType();

					// Update output statistics
					stats.updateRhoOpt( rho );
					stats.updateSolutionStatus( algoStat );

					// Print solution type
					if (options.getPrintLevel() > PrintLevel::NONE)
						MessageHandler::PrintSolution( algoStat );

					return SUCCESSFUL_RETURN;
				} else {
					updatePenalty();

					// Update iterate counters
					updateOuterIter();
					innerIter = -1;
				}
			}

			// Step computation
			ret = solveQPSubproblem( false );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage(ret);
			}

			// Step length computation
			getOptimalStepLength( );

			// Update the total iteration counter
			updateTotalIter();

			// (Failed) termination condition
			if ( totalIter > options.getMaxIterations() )
				return MAX_ITERATIONS_REACHED;

			// Update inner iterate counter
			innerIter++;
		}
	}


	ReturnValue LCQProblem::setConstraints( 	const double* const S1_new, const double* const S2_new,
												const double* const A_new, const double* const lbA_new, const double* const ubA_new )
	{
		if ( nV == 0 || nComp == 0 )
			return LCQPOBJECT_NOT_SETUP;

		if ( A_new == 0 && nC > 0)
			return INVALID_CONSTRAINT_MATRIX;

		// Set up new constraint matrix (A; L; R)
		A = new double[(nC + 2*nComp)*nV];

		for (int i = 0; i < nC*nV; i++)
			A[i] = A_new[i];

		for (int i = 0; i < nComp*nV; i++)
			A[i + nC*nV] = S1_new[i];

		for (int i = 0; i < nComp*nV; i++)
			A[i + nC*nV + nComp*nV] = S2_new[i];

		// Set up new constraint bounds (lbA; 0; 0) & (ubA; INFINITY; INFINITY)
		lbA = new double[nC + 2*nComp];
		ubA = new double[nC + 2*nComp];

		if ( lbA_new != 0 )
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = lbA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = -INFINITY;
		}

		if ( ubA_new != 0 )
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = ubA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = INFINITY;
		}

		for (int i = 0; i < 2*nComp; i++) {
			lbA[i + nC] = 0;
			ubA[i + nC] = INFINITY;
		}

		// Set complementarities
		if ( S1_new == 0 || S2_new == 0 )
			return INVALID_COMPLEMENTARITY_MATRIX;

		S1 = new double[nComp*nV];
		S2 = new double[nComp*nV];

		for (int i = 0; i < nComp*nV; i++) {
			S1[i] = S1_new[i];
			S2[i] = S2_new[i];
		}

		C = new double[nV*nV];
		Utilities::MatrixSymmetrizationProduct(S1, S2, C, nComp, nV);

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::setConstraints(	double* S1_data, int* S1_i, int* S1_p,
											double* S2_data, int* S2_i, int* S2_p,
											double* A_data, int* A_i, int* A_p,
											double* lbA_new, double* ubA_new
											)
	{
		// Create sparse matrices
		S1_sparse = Utilities::copyCSC(nComp, nV, S1_p[nV], S1_data, S1_i, S1_p);
		S2_sparse = Utilities::copyCSC(nComp, nV, S2_p[nV], S2_data, S2_i, S2_p);

		// Get number of elements
		int tmpA_nnx = S1_p[nV] + S2_p[nV];
		tmpA_nnx += A_p != 0 ? A_p[nV] : 0;

		// Data array
		double* tmpA_data = (double*)malloc((size_t)tmpA_nnx*sizeof(double));

		// Row indices
		int* tmpA_i = (int*)malloc((size_t)tmpA_nnx*sizeof(int));

		// Column pointers
		int* tmpA_p = (int*)malloc((size_t)(nV+1)*sizeof(int));

		int index_data = 0;
		tmpA_p[0] = 0;

		// Iterate over columns
		for (int i = 0; i < nV; i++) {
			tmpA_p[i+1] = tmpA_p[i];

			// First handle rows of A
			if (A_p != 0) {
				for (int j = A_p[i]; j < A_p[i+1]; j++) {
					tmpA_data[index_data] = A_data[j];
					tmpA_i[index_data] = A_i[j];
					index_data++;
					tmpA_p[i+1]++;
				}
			}

			// Then rows of S1
			for (int j = S1_p[i]; j < S1_p[i+1]; j++) {
				tmpA_data[index_data] = S1_data[j];
				tmpA_i[index_data] = nC + S1_i[j];
				index_data++;
				tmpA_p[i+1]++;
			}

			// Then rows of S2
			for (int j = S2_p[i]; j < S2_p[i+1]; j++) {
				tmpA_data[index_data] = S2_data[j];
				tmpA_i[index_data] = nC + nComp + S2_i[j];
				index_data++;
				tmpA_p[i+1]++;
			}
		}

		// Create sparse matrix
		A_sparse = Utilities::createCSC(nC + 2*nComp, nV, tmpA_nnx, tmpA_data, tmpA_i, tmpA_p);

		// Set up new constraint bounds (lbA; 0; 0) & (ubA; INFINITY; INFINITY)
		lbA = new double[nC + 2*nComp];
		ubA = new double[nC + 2*nComp];

		if ( lbA_new != 0 )
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = lbA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = -INFINITY;
		}

		if ( ubA_new != 0 )
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = ubA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = INFINITY;
		}

		for (int i = 0; i < 2*nComp; i++) {
			lbA[i + nC] = 0;
			ubA[i + nC] = INFINITY;
		}

		C_sparse = Utilities::MatrixSymmetrizationProduct(S1_data, S1_i, S1_p, S2_data, S2_i, S2_p, nComp, nV);

		if (C_sparse == 0) {
			return FAILED_SYM_COMPLEMENTARITY_MATRIX;
		}

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::setH( double* H_data, int* H_i, int* H_p )
	{
		if (nV <= 0)
			return LCQPOBJECT_NOT_SETUP;

		H_sparse = Utilities::copyCSC(nV, nV, H_p[nV], H_data, H_i, H_p);

		return ReturnValue::SUCCESSFUL_RETURN;
	}


	void LCQProblem::setQk( )
	{
		if (sparseSolver) {
			std::vector<double> Qk_data;
			std::vector<int> Qk_row;
			int* Qk_p = (int*)malloc((size_t)(nV+1)*sizeof(int));
			Qk_p[0] = 0;

			// Iterate over columns
			for (int j = 0; j < nV; j++) {
				Qk_p[j+1] = Qk_p[j];

				int idx_H = H_sparse->p[j];
				int idx_C = C_sparse->p[j];

				while (idx_H < H_sparse->p[j+1] || idx_C < C_sparse->p[j+1]) {
					// Only elements of H left in this column
					if (idx_H < H_sparse->p[j+1] && idx_C >= C_sparse->p[j+1]) {
						Qk_data.push_back(H_sparse->x[idx_H]);
						Qk_row.push_back(H_sparse->i[idx_H]);

						idx_H++;
					}

					// Only elements of C left in this column
					else if (idx_H >= H_sparse->p[j+1] && idx_C < C_sparse->p[j+1]) {
						Qk_data.push_back(rho*C_sparse->x[idx_C]);
						Qk_row.push_back(C_sparse->i[idx_C]);

						Qk_indices_of_C.push_back(Qk_p[j+1]);

						idx_C++;
					}

					// Element of H is in higher row than C
					else if (H_sparse->i[idx_H] < C_sparse->i[idx_C]) {
						Qk_data.push_back(H_sparse->x[idx_H]);
						Qk_row.push_back(H_sparse->i[idx_H]);

						idx_H++;
					}

					// Element of C are in higher row than H
					else if (H_sparse->i[idx_H] > C_sparse->i[idx_C]) {
						Qk_data.push_back(rho*C_sparse->x[idx_C]);
						Qk_row.push_back(C_sparse->i[idx_C]);

						Qk_indices_of_C.push_back(Qk_p[j+1]);

						idx_C++;
					}

					// Add both and update index only once!
					else {
						Qk_data.push_back(H_sparse->x[idx_H] + rho*C_sparse->x[idx_C]);
						Qk_row.push_back(H_sparse->i[idx_H]);

						idx_H++;
						idx_C++;

						Qk_indices_of_C.push_back(Qk_p[j+1]);
					}

					Qk_p[j+1]++;
				}
			}

			int Qk_nnx = (int) Qk_data.size();
			double* Qk_x = (double*)malloc((size_t)Qk_nnx*sizeof(double));
			int* Qk_i =  (int*)malloc((size_t)Qk_nnx*sizeof(int));

			for (size_t i = 0; i < (size_t) Qk_nnx; i++) {
				Qk_x[i] = Qk_data[i];
				Qk_i[i] = Qk_row[i];
			}

			Qk_sparse = Utilities::createCSC(nV, nV, Qk_nnx, Qk_x, Qk_i, Qk_p);

			Qk_data.clear();
			Qk_row.clear();
		} else {
			Utilities::WeightedMatrixAdd(1, H, rho, C, Qk, nV, nV);
		}
	}


	ReturnValue LCQProblem::initializeSolver( )
	{
		ReturnValue ret = SUCCESSFUL_RETURN;
		if (options.getQPSolver() == QPSolver::QPOASES_DENSE) {
			nDuals = nV + nC + 2*nComp;
			boxDualOffset = nV;

			if (sparseSolver) {
				ret = switchToDenseMode( );
				if (ret != SUCCESSFUL_RETURN)
					return ret;
			}

			ret = setLB( lb_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			ret = setUB ( ub_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			Subsolver tmp(nV, nC + 2*nComp, H, A);
			subsolver = tmp;
		} else if (options.getQPSolver() == QPSolver::QPOASES_SPARSE) {
			nDuals = nV + nC + 2*nComp;
			boxDualOffset = nV;

			if (!sparseSolver) {
				ret = switchToSparseMode( );
				if (ret != SUCCESSFUL_RETURN)
					return ret;
			}

			Subsolver tmp(nV, nC + 2*nComp, H_sparse, A_sparse, options.getQPSolver());

			ret = setLB( lb_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			ret = setUB( ub_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;


			subsolver = tmp;
		} else if (options.getQPSolver() == QPSolver::OSQP_SPARSE) {
			if (lb != 0 || ub != 0) {
				return INVALID_OSQP_BOX_CONSTRAINTS;
			}

			nDuals = nC + 2*nComp;
			boxDualOffset = 0;

			if (!sparseSolver) {
				ret = switchToSparseMode( );
				if (ret != SUCCESSFUL_RETURN)
					return ret;
			}

			if (lb_tmp != 0 || ub_tmp != 0) {
				return ReturnValue::INVALID_OSQP_BOX_CONSTRAINTS;
			}

			Subsolver tmp(nV, nDuals, H_sparse, A_sparse, options.getQPSolver());
			subsolver = tmp;
		} else {
			return ReturnValue::NOT_YET_IMPLEMENTED;
		}

		// Initialize variables and counters
		alphak = 1;
		rho = options.getInitialPenaltyParameter( );
		outerIter = 0;
		innerIter = 0;
		totalIter = 0;
		algoStat = AlgorithmStatus::PROBLEM_NOT_SOLVED;

		// Initialize subproblem solver
		subsolver.setPrintLevel( options.getPrintLevel() );

		// Reset output statistics
		stats.reset();

		// Set seed
		srand( (unsigned int)time( NULL ) );

		// Clean up temporary box constraints
		if (lb_tmp != 0) {
			delete[] lb_tmp;
			lb_tmp = NULL;
		}

		if (ub_tmp != 0) {
			delete[] ub_tmp;
			ub_tmp = NULL;
		}

		return ret;
	}


	ReturnValue LCQProblem::switchToSparseMode( )
	{
		H_sparse = Utilities::dns_to_csc(H, nV, nV);
		A_sparse = Utilities::dns_to_csc(A, nC + 2*nComp, nV);

		S1_sparse = Utilities::dns_to_csc(S1, nComp, nV);
		S2_sparse = Utilities::dns_to_csc(S2, nComp, nV);
		C_sparse = Utilities::dns_to_csc(C, nV, nV);

		// Make sure that all sparse matrices are not null pointer
		if (H_sparse == 0 || A_sparse == 0 || S1_sparse == 0 || S2_sparse == 0 || C_sparse == 0) {
			return FAILED_SWITCH_TO_SPARSE;
		}

		// Clean up dense data (only if succeeded)
		delete[] H; H = NULL;
		delete[] A; A = NULL;
		delete[] S1; S1 = NULL;
		delete[] S2; S2 = NULL;
		delete[] C; C = NULL;

		// Toggle sparsity flag
		sparseSolver = true;

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::switchToDenseMode( )
	{
		H = Utilities::csc_to_dns(H_sparse);
		A = Utilities::csc_to_dns(A_sparse);

		S1 = Utilities::csc_to_dns(S1_sparse);
		S2 = Utilities::csc_to_dns(S2_sparse);
		C = Utilities::csc_to_dns(C_sparse);

		// Make sure that all sparse matrices are not null pointer
		if (H == 0 || A == 0 || S1 == 0 || S2 == 0 || C == 0) {
			return FAILED_SWITCH_TO_DENSE;
		}

		// Clean up sparse data (only if succeeded)
		Utilities::ClearSparseMat(&C_sparse);
		Utilities::ClearSparseMat(&A_sparse);
		Utilities::ClearSparseMat(&H_sparse);
		Utilities::ClearSparseMat(&S1_sparse);
		Utilities::ClearSparseMat(&S2_sparse);

		// Toggle sparsity flag
		sparseSolver = false;

		return SUCCESSFUL_RETURN;
	}


	void LCQProblem::updateLinearization()
	{
		if (sparseSolver) {
			Utilities::AffineLinearTransformation(rho, C_sparse, xk, g, gk, nV);
		} else {
			Utilities::AffineLinearTransformation(rho, C, xk, g, gk, nV, nV);
		}
	}


	ReturnValue LCQProblem::solveQPSubproblem(bool initialSolve)
	{
		// First solve convex subproblem
		ReturnValue ret = subsolver.solve( initialSolve, qpIterk, gk, lbA, ubA, xk, yk, lb, ub );

		// Update stats
		stats.updateSubproblemIter(qpIterk);

		// If no initial guess was passed, then need to allocate memory
		if (xk == 0) {
			xk = new double[nV];
		}

		if (yk == 0) {
			yk = new double[nDuals];
		}

		// Return on error
		if (ret != SUCCESSFUL_RETURN)
			return ret;

		// Update xnew, yk
		subsolver.getSolution(xnew, yk);

		// Update yk_A
		for (int i = 0; i < nC + 2*nComp; i++)
			yk_A[i] = yk[boxDualOffset + i];

		// Update pk
		Utilities::WeightedVectorAdd(1, xnew, -1, xk, pk, nV);

		return SUCCESSFUL_RETURN;
	}


	bool LCQProblem::stationarityCheck( ) {
		return Utilities::MaxAbs(statk, nV) < options.getStationarityTolerance();
	}


	bool LCQProblem::complementarityCheck( ) {
		if (sparseSolver) {
			return Utilities::QuadraticFormProduct(C_sparse, xk, nV) < 2*options.getComplementarityTolerance();
		} else {
			return Utilities::QuadraticFormProduct(C, xk, nV) < 2*options.getComplementarityTolerance();
		}

	}


	void LCQProblem::updatePenalty( ) {
		rho *= options.getPenaltyUpdateFactor();
	}


	void LCQProblem::getOptimalStepLength( ) {

		double qk;

		if (sparseSolver) {
			qk = Utilities::QuadraticFormProduct(Qk_sparse, pk, nV);
			Utilities::AffineLinearTransformation(1, Qk_sparse, xk, g, lk_tmp, nV);
		} else {
			qk = Utilities::QuadraticFormProduct(Qk, pk, nV);
			Utilities::AffineLinearTransformation(1, Qk, xk, g, lk_tmp, nV, nV);
		}

		double lk = Utilities::DotProduct(pk, lk_tmp, nV);

		alphak = 1;

		// Convex Descent Case
		if (qk > 0 && lk < 0) {
			alphak = std::min(-lk/qk, 1.0);
		}
	}


	void LCQProblem::updateStep( ) {
		// Update penalty on rejected step
		// Currently doesn't exist
		if (alphak <= 0) {
			updatePenalty( );

			// Update iterate counters
			updateOuterIter();
			innerIter = 0;
		} else {
			// xk = xk + alphak*pk
			Utilities::WeightedVectorAdd(1, xk, alphak, pk, xk, nV);

			// gk = new linearization + g
			updateLinearization();

			// Add some +/- EPS to each coordinate
			perturbGradient();
		}

		// Update Qk = H + rho*C (only required on first inner iteration)
		if (innerIter == 0 && outerIter == 0) {
			setQk();
		} else if (innerIter == 0) {
			updateQk();
		}

		// stat = Qk*xk + g - A'*yk_A - yk_x
		// 1) Objective contribution: Qk*xk + g
		if (sparseSolver) {
			Utilities::AffineLinearTransformation(1, Qk_sparse, xk, g, statk, nV);
		} else {
			Utilities::AffineLinearTransformation(1, Qk, xk, g, statk, nV, nV);
		}

		// 2) Constraint contribution: A'*yk
		// Utilities::TransponsedMatrixMultiplication(A, yk_A, constr_statk, nC + 2*nComp, nV, 1);
		if (sparseSolver) {
			Utilities::TransponsedMatrixMultiplication(A_sparse, yk_A, constr_statk, nC + 2*nComp, nV);
		} else {
			Utilities::TransponsedMatrixMultiplication(A, yk_A, constr_statk, nC + 2*nComp, nV, 1);
		}

		Utilities::WeightedVectorAdd(1, statk, -1, constr_statk, statk, nV);

		// 3) Box constraint contribution
		if (lb != 0 || ub != 0) {
			for (int i = 0; i < nV; i++)
				box_statk[i] = yk[i];

			Utilities::WeightedVectorAdd(1, statk, -1, box_statk, statk, nV);
		}
	}


	void LCQProblem::updateQk( ) {
		// Smart update in sparse case
		if (sparseSolver) {
			double factor = rho*(1 - 1.0/options.getPenaltyUpdateFactor());
			for (size_t j = 0; j < Qk_indices_of_C.size(); j++) {
				Qk_sparse->x[Qk_indices_of_C[j]] += factor*C_sparse->x[j];
			}
		} else {
			Utilities::WeightedMatrixAdd(1, H, rho, C, Qk, nV, nV);
		}

	}


	void LCQProblem::updateOuterIter( ) {
		outerIter++;
		stats.updateIterOuter(1);
	}


	void LCQProblem::updateTotalIter( ) {
		totalIter++;
		stats.updateIterTotal(1);
	}


	void LCQProblem::perturbGradient( ) {

		int randNum;
		for (int i = 0; i < nV; i++) {
			// Random number -1, 0, 1
			randNum = (rand() % 3) - 1;

			gk[i] += randNum*Utilities::EPS;
		}
	}


	void LCQProblem::perturbStep( ) {

		int randNum;
		for (int i = 0; i < nV; i++) {
			// Random number -1, 0, 1
			randNum = (rand() % 3) - 1;

			xk[i] += randNum*Utilities::EPS;
		}
	}


	void LCQProblem::transformDuals( ) {

		double* tmp = new double[nComp];

		// y_S1 = y - rho*S2*xk
		if (sparseSolver) {
			Utilities::MatrixMultiplication(S2_sparse, xk, tmp);
		} else {
			Utilities::MatrixMultiplication(S2, xk, tmp, nComp, nV, 1);
		}

		for (int i = 0; i < nComp; i++) {
			yk[boxDualOffset + nC + i] = yk[boxDualOffset + nC + i] - rho*tmp[i];
		}

		// y_S2 = y - rho*S1*xk
		if (sparseSolver) {
			Utilities::MatrixMultiplication(S1_sparse, xk, tmp);
		} else {
			Utilities::MatrixMultiplication(S1, xk, tmp, nComp, nV, 1);
		}

		for (int i = 0; i < nComp; i++) {
			yk[boxDualOffset + nC + nComp + i] = yk[boxDualOffset + nC + nComp + i] - rho*tmp[i];
		}

		// clear memory
		delete[] tmp;
	}


	void LCQProblem::determineStationarityType( ) {

		std::vector<int> weakComp = getWeakComplementarities( );

		bool s_stationary = true;
		bool m_stationary = true;

		for (int i : weakComp) {
			double dualProd = yk_A[nC + i]*yk_A[nC + nComp + i];
			double dualMin = std::min(yk_A[nC + i], yk_A[nC + nComp + i]);

			// Check failure of s-stationarity
			if (dualMin < 0)
				s_stationary = false;

			// Check failure of m-/c-stationarity
			if (std::abs(dualProd) >= options.getComplementarityTolerance() && dualMin <= 0) {

				// Check failure of c-stationarity
				if (dualProd <= options.getComplementarityTolerance()) {
					algoStat = AlgorithmStatus::W_STATIONARY_SOLUTION;
					return;
				}

				m_stationary = false;
			}
		}

		if (s_stationary) {
			algoStat = AlgorithmStatus::S_STATIONARY_SOLUTION;
			return;
		}


		if (m_stationary) {
			algoStat = AlgorithmStatus::M_STATIONARY_SOLUTION;
			return;
		}

		algoStat = AlgorithmStatus::C_STATIONARY_SOLUTION;
		return;
	}


	std::vector<int> LCQProblem::getWeakComplementarities( )
	{
		double* S1x = new double[nComp];
		double* S2x = new double[nComp];

		if (sparseSolver) {
			Utilities::MatrixMultiplication(S1_sparse, xk, S1x);
			Utilities::MatrixMultiplication(S2_sparse, xk, S2x);
		} else {
			Utilities::MatrixMultiplication(S1, xk, S1x, nComp, nV, 1);
			Utilities::MatrixMultiplication(S2, xk, S2x, nComp, nV, 1);
		}

		std::vector<int> indices;

		for (int i = 0; i < nComp; i++) {
			if (S1x[i] <= options.getComplementarityTolerance())
				if (S2x[i] <= options.getComplementarityTolerance())
					indices.push_back(i);
		}

		// Free memory
		delete[] S1x;
		delete[] S2x;

		return indices;
	}


	AlgorithmStatus LCQProblem::getPrimalSolution( double* const xOpt ) const
	{
		if (algoStat != AlgorithmStatus::PROBLEM_NOT_SOLVED) {
			for (int i = 0; i < nV; i++)
				xOpt[i] = xk[i];
		}

		return algoStat;
	}


	AlgorithmStatus LCQProblem::getDualSolution( double* const yOpt ) const
	{
		if (algoStat != AlgorithmStatus::PROBLEM_NOT_SOLVED) {
			for (int i = 0; i < nDuals; i++)
				yOpt[i] = yk[i];
		}

		return algoStat;
	}


	int LCQProblem::getNumerOfDuals( ) const
	{
		return nDuals;
	}


	void LCQProblem::getOutputStatistics( OutputStatistics& _stats) const
	{
		_stats = stats;
	}


	/*
	 *	 p r i n t I t e r a t i o n
	 */
	void LCQProblem::printIteration( )
	{
		if (options.getPrintLevel() == PrintLevel::NONE)
			return;

		if (options.getPrintLevel() == PrintLevel::OUTER_LOOP_ITERATES && innerIter % 10 > 0)
			return;

		// Print header every 10 iters
		bool headerInner = (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES && innerIter % 10 == 0);
		bool headerOuter = (options.getPrintLevel() == PrintLevel::OUTER_LOOP_ITERATES && outerIter % 10 == 0);

		if (headerInner || headerOuter)
			printHeader();

		const char* sep = " | ";

		// Print outer iterate
		printf("%6d", outerIter);

		// Print innter iterate
		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES)
			printf("%s%*d", sep, 6, innerIter);

		// Print stationarity violation
		double tmpdbl = Utilities::MaxAbs(statk, nV);
		printf("%s%10.3g", sep, tmpdbl);

		// Print complementarity violation
		if (sparseSolver) {
			tmpdbl = Utilities::QuadraticFormProduct(C_sparse, xk, nV)/2.0;
		} else {
			tmpdbl = Utilities::QuadraticFormProduct(C, xk, nV)/2.0;
		}

		printf("%s%10.3g", sep, tmpdbl);

		// Print current penalty parameter
		printf("%s%10.3g", sep, rho);

		// Print infinity norm of computed full step
		tmpdbl = Utilities::MaxAbs(pk, nV);
		printf("%s%10.3g", sep, tmpdbl);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES) {
			// Print optimal step length
			printf("%s%10.3g", sep, alphak);

			// Print number of qpOASES iterations
			printf("%s%6d", sep, qpIterk);
		}

		// Print new line
		printf(" \n");
	}


	void LCQProblem::printHeader()
	{
		printLine();

		const char* sep = " | ";
		const char* outer = " outer";
		const char* inner = " inner";
		const char* stat = "  station ";
		const char* comp = "  complem ";
		const char* pen = "    rho   ";
		const char* np = "  norm p  ";
		const char* sl = "   alpha  ";
		const char* subIt = "sub it";

		printf("%s",outer);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES)
			printf("%s%s", sep, inner);

		printf("%s%s", sep, stat);
		printf("%s%s", sep, comp);
		printf("%s%s", sep, pen);
		printf("%s%s", sep, np);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES) {
			printf("%s%s%s%s", sep, sl, sep, subIt);
		}

		printf(" \n");

		printLine();
	}

	/*
	 * 	 p r i n t L i n e
	 */
	void LCQProblem::printLine()
	{
		const char* iSep = "------";
		const char* dSep = "----------";
		const char* node = "-+-";

		if (innerIter == 0 && outerIter == 0) {
			// Print new line at very beginning as there is no guarantee
			// that whatever was printed before ended with newline.
			printf("\n");
		}

		printf("%s", iSep);

		// Print innter iterate
		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES)
			printf("%s%s", node, iSep);

		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES) {
			printf("%s%s%s%s", node, dSep, node, iSep);
		}

		printf("-\n");
	}


	/// Clear allocated memory
	void LCQProblem::clear( )
	{
		if (H != 0) {
			delete[] H;
			H = NULL;
		}

		if (g != 0) {
			delete[] g;
			g = NULL;
		}

		if (lb != 0) {
			delete[] lb;
			lb = NULL;
		}

		if (ub != 0) {
			delete[] ub;
			ub = NULL;
		}

		if (lb_tmp != 0) {
			delete[] lb_tmp;
			lb_tmp = NULL;
		}

		if (ub_tmp != 0) {
			delete[] ub_tmp;
			ub_tmp = NULL;
		}

		if (A != 0) {
			delete[] A;
			A = NULL;
		}

		if (lbA != 0) {
			delete[] lbA;
			lbA = NULL;
		}

		if (ubA != 0) {
			delete[] ubA;
			ubA = NULL;
		}

		if (S1 != 0) {
			delete[] S1;
			S1 = NULL;
		}

		if (S2 != 0) {
			delete[] S2;
			S2 = NULL;
		}

		if (C != 0) {
			delete[] C;
			C = NULL;
		}

		if (gk != 0) {
			delete[] gk;
			gk = NULL;
		}

		if (xk != 0) {
			delete[] xk;
			xk = NULL;
		}

		if (yk != 0) {
			delete[] yk;
			yk = NULL;
		}

		if (yk_A != 0) {
			delete[] yk_A;
			yk_A = NULL;
		}

		if (xnew != 0) {
			delete[] xnew;
			xnew = NULL;
		}

		if (pk != 0) {
			delete[] pk;
			pk = NULL;
		}

		if (Qk != 0) {
			delete[] Qk;
			Qk = NULL;
		}

		if (statk != 0) {
			delete[] statk;
			statk = NULL;
		}

		if (constr_statk != 0) {
			delete[] constr_statk;
			constr_statk = NULL;
		}

		if (box_statk != 0) {
			delete[] box_statk;
			box_statk = NULL;
		}

		if (lk_tmp != 0) {
			delete[] lk_tmp;
			lk_tmp = NULL;
		}

		Utilities::ClearSparseMat(&C_sparse);
		Utilities::ClearSparseMat(&A_sparse);
		Utilities::ClearSparseMat(&H_sparse);
		Utilities::ClearSparseMat(&Qk_sparse);
		Utilities::ClearSparseMat(&S1_sparse);
		Utilities::ClearSparseMat(&S2_sparse);
	}
}





/*
 *	end of file
 */