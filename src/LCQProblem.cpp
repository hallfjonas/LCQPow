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


#include "LCQProblem.hpp"
#include "Utilities.hpp"
#include "SubsolverQPOASES.hpp"
#include "SubsolverOSQP.hpp"
#include "PlotManager.hpp"

#include <iostream>
#include <string>
#include <math.h>
#include <qpOASES.hpp>

using qpOASES::QProblem;

namespace lcqpOASES {

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


	returnValue LCQProblem::solve(	const double* const _H, const double* const _g,
									const double* const _lb, const double* const _ub,
									const double* const _S1, const double* const _S2,
									const double* const _A, const double* const _lbA, const double* const _ubA,
									const double* const _x0, const double* const _y0
									)
	{
		returnValue ret;

		ret = setH( _H );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setLB( _lb );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setUB( _ub );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setConstraints( _S1, _S2, _A, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		// Number of duals in qpOASES case (box constraints + linear constraints + complementarity bounds):
		nDuals = nV + nC + 2*nComp;
		boxDualOffset = nV;

		// USe qpOASES in dense formulations
		Subsolver tmp( nV, nC + 2*nComp, H, A );
		subsolver = tmp;

		return MessageHandler::PrintMessage( runSolver( ) );
	}


	returnValue LCQProblem::solve(	const char* const H_file, const char* const g_file,
									const char* const lb_file, const char* const ub_file,
									const char* const S1_file, const char* const S2_file,
									const char* const A_file, const char* const lbA_file, const char* const ubA_file,
									const char* const x0_file, const char* const y0_file
									)
	{
		returnValue ret;

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

		double* _lb = new double[nV];
		ret = Utilities::readFromFile( _lb, nV, lb_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _lb;
			return MessageHandler::PrintMessage( ret );
		}

		double* _ub = new double[nV];
		ret = Utilities::readFromFile( _ub, nV, ub_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _ub;
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

		double* _x0 = NULL;
		if (x0_file != 0) {
			_x0 = new double[nC + 2*nComp];
			ret = Utilities::readFromFile( _x0, nC + 2*nComp, x0_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _S1;
				return MessageHandler::PrintMessage( ret );
			}
		}

		double* _y0 = NULL;
		if (y0_file != 0) {
			_y0 = new double[nC + 2*nComp];
			ret = Utilities::readFromFile( _y0, nC + 2*nComp, y0_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _S2;
				return MessageHandler::PrintMessage( ret );
			}
		}

		// Fill vaues
		ret = setH( _H );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setLB( _lb );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setUB( _ub );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setConstraints( _S1, _S2, _A, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		// Number of duals in qpOASES case (box constraints + linear constraints + complementarity bounds):
		nDuals = nV + nC + 2*nComp;
		boxDualOffset = nV;

		// USe qpOASES in dense formulations
		Subsolver tmp( nV, nC + 2*nComp, H, A );
		subsolver = tmp;

		// Call solver
		return MessageHandler::PrintMessage( runSolver( ) );
	}


	returnValue LCQProblem::solve(	double* _H_data, int _H_nnx, int* _H_i, int* _H_p, double* _g,
									double* _S1_data, int _S1_nnx, int* _S1_i, int* _S1_p,
									double* _S2_data, int _S2_nnx, int* _S2_i, int* _S2_p,
									double* _A_data, int _A_nnx, int* _A_i, int* _A_p,
									double* _lbA, double* _ubA,	double* _x0, double* _y0
									)
	{
		returnValue ret;

		ret = setH( _H_data, _H_nnx, _H_i, _H_p );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setConstraints( _S1_data, _S1_nnx, _S1_i, _S1_p, _S2_data, _S2_nnx, _S2_i, _S2_p, _A_data, _A_nnx, _A_i, _A_p, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret );

		// Make sure that box constraints are null pointers in OSQP case
		lb = NULL;
		ub = NULL;

		// Number of duals in OSQP case:
		nDuals = nC + 2*nComp;
		boxDualOffset = 0;

		// Use OSQP in sparse formulations
		Subsolver tmp(nV, nC + 2*nComp, H_sparse, A_sparse, g, lbA, ubA);
		subsolver = tmp;

		return MessageHandler::PrintMessage( runSolver( ) );

	}


	returnValue LCQProblem::setConstraints( 	const double* const S1_new, const double* const S2_new,
												const double* const A_new, const double* const lbA_new, const double* const ubA_new )
	{
		if ( nV == 0 || nComp == 0 )
			return LCQPOBJECT_NOT_SETUP;

		if ( A_new == 0 && nC > 0)
			return INVALID_ARGUMENT;

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
			return INVALID_ARGUMENT;

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


	returnValue LCQProblem::setConstraints(	double* S1_data, int S1_nnx, int* S1_i, int* S1_p,
											double* S2_data, int S2_nnx, int* S2_i, int* S2_p,
											double* A_data, int A_nnx, int* A_i, int* A_p,
											double* lbA_new, double* ubA_new
											)
	{

		int tmpA_nnx = A_nnx + S1_nnx + S2_nnx;
		tmpA_data = new double[tmpA_nnx];
		tmpA_i = new int[tmpA_nnx];
		tmpA_p = new int[nV+1];

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
		A_sparse = csc_matrix(nC + 2*nComp, nV, tmpA_nnx, tmpA_data, tmpA_i, tmpA_p);

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


		// Create dense matrices
		A = new double[(nC + 2*nComp)*nV]();
		lcqpOASES::returnValue ret = Utilities::csc_to_dns(A_sparse, A, nC + 2*nComp, nV);
		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage(ret);

		// For now store S1, S2, and C in dense format.
		// If we can manage to adapt all operations required for C (specifically S1'*S2 + S2'*S1 = C)
		// we should instantly switch to sparse format!
		S1 = new double[nComp*nV]();
		for (int j = 0; j < nV; j++) {
			for (int i = 0; i < S1_p[j+1] - S1_p[j]; i++) {
				S1[(S1_p[j]+i)*nV + j] = S1_data[S1_p[j]+i];
			}
		}

		S2 = new double[nComp*nV]();
		for (int j = 0; j < nV; j++) {
			for (int i = 0; i < S2_p[j+1] - S2_p[j]; i++) {
				S2[(S2_p[j]+i)*nV + j] = S2_data[S2_p[j]+i];
			}
		}

		C = new double[nV*nV];
		Utilities::MatrixSymmetrizationProduct(S1, S2, C, nComp, nV);

		return SUCCESSFUL_RETURN;
	}


	returnValue LCQProblem::setH( double* H_data, int H_nnx, int* H_i, int* H_p )
	{
		if (nV <= 0)
			return LCQPOBJECT_NOT_SETUP;

		H_sparse = csc_matrix(nV, nV, H_nnx, H_data, H_i, H_p);

		H = new double[nV*nV]();
		Utilities::csc_to_dns(H_sparse, H, nV, nV);

		return returnValue::SUCCESSFUL_RETURN;
	}


	returnValue LCQProblem::runSolver( )
	{
		// Create a plot manager instance
		PlotManager plotter(nV, nC, nComp, LCQPNAME::IVOCP);

		// Initialize variables
		initializeSolver();

		// Initialization strategy
		if (options.solveZeroPenaltyFirst) {
			memcpy(gk, g, (long unsigned int)nV*sizeof(double));

			if (solveQPSubproblem( true ) != SUCCESSFUL_RETURN) {
				return INITIAL_SUBPROBLEM_FAILED;
			}

			Utilities::AffineLinearTransformation(rho, C, xk, g, gk, nV, nV);

			if (solveQPSubproblem( false ) != SUCCESSFUL_RETURN) {
				return SUBPROBLEM_SOLVER_ERROR;
			}

		} else {
			Utilities::AffineLinearTransformation(rho, C, xk, g, gk, nV, nV);

			if (solveQPSubproblem( true ) != SUCCESSFUL_RETURN) {
				return INITIAL_SUBPROBLEM_FAILED;
			}
		}

		// Outer and inner loop in one
		while ( true ) {

			// Update xk, gk, Qk, stationarity
			updateStep( );

			// Create debugging plots
			// plotter.CreateIVOCPPlots(xk, lb, ub);

			// Print iteration
			printIteration( );

			// Terminate, update pen, or continue inner loop
			if (stationarityCheck()) {
				if (complementarityCheck()) {
					// Switch from penalized to LCQP duals
					transformDuals();

					// Determine C-,M-,S-Stationarity
					determineStationarityType();

					// Print solution type
					if (options.printLvl > printLevel::NONE)
						MessageHandler::PrintSolution( algoStat );

					return SUCCESSFUL_RETURN;
				} else {
					updatePenalty();

					// Update iterate counters
					outerIter++;
					innerIter = -1;
				}
			}

			// Step computation
			if (solveQPSubproblem( false ) != SUCCESSFUL_RETURN) {
				return SUBPROBLEM_SOLVER_ERROR;
			}

			// Step length computation
			getOptimalStepLength( );

			if ( outerIter > options.maxOuterIterations )
				return MAX_OUTER_ITERATIONS_REACHED;

			if ( innerIter > options.maxInnerIterations )
				return MAX_INNER_ITERATIONS_REACHED;

			// Update inner iterate counter
			innerIter++;
		}
	}


	void LCQProblem::initializeSolver( )
	{
		// Initialize variables and counters
		alphak = 1;
		rho = options.initialComplementarityPenalty;
		outerIter = 0;
		innerIter = 0;
		algoStat = algorithmStatus::PROBLEM_NOT_SOLVED;

		// Initialize subproblem solver
		subsolver.setPrintLevel( options.printLvl );
	}


	returnValue LCQProblem::solveQPSubproblem(bool initialSolve)
	{
		// First solve convex subproblem
		returnValue ret = subsolver.solve( initialSolve, qpIterk, gk, lbA, ubA, xk, yk, lb, ub );

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
		return Utilities::MaxAbs(statk, nV) < options.stationarityTolerance;
	}


	bool LCQProblem::complementarityCheck( ) {
		return Utilities::QuadraticFormProduct(C, xk, nV) < 2*options.complementarityTolerance;
	}


	void LCQProblem::updatePenalty( ) {
		rho *= options.complementarityPenaltyUpdate;
	}


	void LCQProblem::getOptimalStepLength( ) {

		double qk = Utilities::QuadraticFormProduct(Qk, pk, nV);

		Utilities::AffineLinearTransformation(1, Qk, xk, g, lk_tmp, nV, nV);

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
			outerIter++;
			innerIter = 0;
		} else {
			// xk = xk + alphak*pk
			Utilities::WeightedVectorAdd(1, xk, alphak, pk, xk, nV);

			// gk = new linearization + g
			Utilities::AffineLinearTransformation(rho, C, xk, g, gk, nV, nV);

			// Add some +/- EPS to each coordinate
			perturbGradient();
			// perturbStep();
		}

		// Update Qk = H + rho*C (only required on first inner iteration)
		if (innerIter == 0)
			Utilities::WeightedMatrixAdd(1, H, rho, C, Qk, nV, nV);

		// stat = Qk*xk + g - A'*yk_A - yk_x
		// 1) Objective contribution: Qk*xk + g
		Utilities::AffineLinearTransformation(1, Qk, xk, g, statk, nV, nV);

		// 2) Constraint contribution: A'*yk
		Utilities::TransponsedMatrixMultiplication(A, yk_A, constr_statk, nC + 2*nComp, nV, 1);
		Utilities::WeightedVectorAdd(1, statk, -1, constr_statk, statk, nV);

		// 3) Box constraint contribution
		if (lb != 0 || ub != 0) {
			for (int i = 0; i < nV; i++)
				box_statk[i] = yk[i];

			Utilities::WeightedVectorAdd(1, statk, -1, box_statk, statk, nV);
		}
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
		Utilities::MatrixMultiplication(S2, xk, tmp, nComp, nV, 1);
		for (int i = 0; i < nComp; i++) {
			yk[boxDualOffset + nC + i] = yk[boxDualOffset + nC + i] - rho*tmp[i];
		}

		// y_S2 = y - rho*S1*xk
		Utilities::MatrixMultiplication(S1, xk, tmp, nComp, nV, 1);
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
			if (std::abs(dualProd) >= options.complementarityTolerance && dualMin <= 0) {

				// Check failure of c-stationarity
				if (dualProd <= options.complementarityTolerance) {
					algoStat = algorithmStatus::W_STATIONARY_SOLUTION;
					return;
				}

				m_stationary = false;
			}
		}

		if (s_stationary) {
			algoStat = algorithmStatus::S_STATIONARY_SOLUTION;
			return;
		}


		if (m_stationary) {
			algoStat = algorithmStatus::M_STATIONARY_SOLUTION;
			return;
		}

		algoStat = algorithmStatus::C_STATIONARY_SOLUTION;
		return;
	}


	std::vector<int> LCQProblem::getWeakComplementarities( )
	{
		double* S1x = new double[nComp];
		double* S2x = new double[nComp];

		Utilities::MatrixMultiplication(S1, xk, S1x, nComp, nV, 1);
		Utilities::MatrixMultiplication(S2, xk, S2x, nComp, nV, 1);

		std::vector<int> indices;

		for (int i = 0; i < nComp; i++) {
			if (S1x[i] <= options.complementarityTolerance)
				if (S2x[i] <= options.complementarityTolerance)
					indices.push_back(i);
		}

		// Free memory
		delete[] S1x;
		delete[] S2x;

		return indices;
	}


	algorithmStatus LCQProblem::getPrimalSolution( double* const xOpt ) const
	{
		if (algoStat != algorithmStatus::PROBLEM_NOT_SOLVED) {
			for (int i = 0; i < nV; i++)
				xOpt[i] = xk[i];
		}

		return algoStat;
	}


	algorithmStatus LCQProblem::getDualSolution( double* const yOpt ) const
	{
		if (algoStat != algorithmStatus::PROBLEM_NOT_SOLVED) {
			for (int i = 0; i < nDuals; i++)
				yOpt[i] = yk[i];
		}

		return algoStat;
	}


	/*
	 *	 p r i n t I t e r a t i o n
	 */
	void LCQProblem::printIteration( )
	{
		if (options.printLvl == printLevel::NONE)
			return;

		// Print header every 10 iters
		bool headerInner = (options.printLvl >= printLevel::INNER_LOOP_ITERATES && innerIter % 10 == 0);
		bool headerOuter = (options.printLvl == printLevel::OUTER_LOOP_ITERATES && outerIter % 10 == 0);

		if (headerInner || headerOuter)
			printHeader();

		const char* sep = " | ";

		// Print outer iterate
		printf("%6d", outerIter);

		// Print innter iterate
		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			printf("%s%*d", sep, 6, innerIter);

		// Print stationarity violation
		double tmpdbl = Utilities::MaxAbs(statk, nV);
		printf("%s%10.3g", sep, tmpdbl);

		// Print complementarity violation
		tmpdbl = Utilities::QuadraticFormProduct(C, xk, nV)/2.0;
		printf("%s%10.3g", sep, tmpdbl);

		// Print current penalty parameter
		printf("%s%10.3g", sep, rho);

		// Print infinity norm of computed full step
		tmpdbl = Utilities::MaxAbs(pk, nV);
		printf("%s%10.3g", sep, tmpdbl);

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
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

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			printf("%s%s", sep, inner);

		printf("%s%s", sep, stat);
		printf("%s%s", sep, comp);
		printf("%s%s", sep, pen);
		printf("%s%s", sep, np);

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
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

		printf("%s", iSep);

		// Print innter iterate
		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			printf("%s%s", node, iSep);

		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
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

		if (tmpA_data != 0) {
			delete[] tmpA_data;
			tmpA_data = NULL;
		}

		if (tmpA_p != 0) {
			delete[] tmpA_p;
			tmpA_p = NULL;
		}

		if (tmpA_i != 0) {
			delete[] tmpA_i;
			tmpA_i = NULL;
		}

		if (H_sparse != 0) {
			c_free(H_sparse);
			H_sparse = NULL;
		}

		if (A_sparse != 0) {
			c_free(A_sparse);
			A_sparse = NULL;
		}
	}
}





/*
 *	end of file
 */