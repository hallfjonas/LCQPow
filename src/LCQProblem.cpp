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
#include <iostream>
#include <string>
#include <math.h>
#include <qpOASES.hpp>

using qpOASES::QProblem;

namespace lcqpOASES {

	/*****************************************************************************
	 *  P U B L I C                                                              *
	 *****************************************************************************/

	/*
	*	Q P r o b l e m
	*/
	LCQProblem::LCQProblem( int _nV, int _nC, int _nComp )
	{
		/* consistency checks */
		if ( _nV <= 0 || _nComp <= 0)
		{
			throw( INVALID_ARGUMENT );
		}

		if ( _nC < 0 )
		{
			_nC = 0;
			throw( INVALID_ARGUMENT );
		}

		nV = _nV;
		nC = _nC;
		nComp = _nComp;

		QProblem lp(nV, nC, qpOASES::HessianType::HST_POSDEF);
		qp = lp;
	}


	/*
	*	s o l v e
	*/
	returnValue LCQProblem::solve(	const double* const _H, const double* const _g, const double* const _A,
									const double* const _lb, const double* const _ub,
									const double* const _lbA, const double* const _ubA,
									const double* const _S1, const double* const _S2, const double* const _x0, const double* const _y0
									)
	{
		returnValue ret = setupLCQPdata(_H, _g, _A, _lb, _ub, _lbA, _ubA, _S1, _S2, _x0, _y0);

		if (ret != SUCCESSFUL_RETURN)
			return ret;

		return runSolver( );
	}


	/*
	*	i n i t
	*/
	returnValue LCQProblem::solve(	const char* const H_file, const char* const g_file, const char* const A_file,
									const char* const lb_file, const char* const ub_file,
									const char* const lbA_file, const char* const ubA_file,
									const char* const S1_file, const char* const S2_file,
									const char* const x0_file, const char* const y0_file
									)
	{
		// Internally set up data
		returnValue ret = setupLCQPdata(H_file, g_file, A_file, lb_file, ub_file, lbA_file, ubA_file, S1_file, S2_file, x0_file, y0_file);

		if (ret != SUCCESSFUL_RETURN)
			return ret;

		return runSolver( );	
	}


	/*
	*	s e t u p Q P d a t a
	*/
	returnValue LCQProblem::setupLCQPdata(	const double* const _H, const double* const _g, const double* const _A,
											const double* const _lb, const double* const _ub,
											const double* const _lbA, const double* const _ubA,
											const double* const _S1, const double* const _S2,
											const double* const _x0, const double* const _y0
											)
	{
		returnValue ret;

		ret = setH( _H );

		if (ret != SUCCESSFUL_RETURN)
			return ret;	

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return ret;	

		ret = setLB( _lb );

		if (ret != SUCCESSFUL_RETURN)
			return ret;	

		ret = setUB( _ub );

		if (ret != SUCCESSFUL_RETURN)
			return ret;	

		ret = setConstraints( _A, _S1, _S2, _lbA, _ubA );	
		
		if (ret != SUCCESSFUL_RETURN)
			return ret;	
		
		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return ret;	

		return SUCCESSFUL_RETURN;
	}

	/*
	*	s e t u p Q P d a t a F r o m F i l e
	*/
	returnValue LCQProblem::setupLCQPdata(	const char* const H_file, const char* const g_file, const char* const A_file,
											const char* const lb_file, const char* const ub_file,
											const char* const lbA_file, const char* const ubA_file,
											const char* const S1_file, const char* const S2_file,
											const char* const x0_file, const char* const y0_file
											)
	{
		returnValue returnvalue;
		
		double* _H = new double[nV*nV];
		returnvalue = Utilities::readFromFile( _H, nV*nV, H_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _H;
			return returnvalue;
		}
			

		double* _g = new double[nV];
		returnvalue = Utilities::readFromFile( _g, nV, g_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _g;
			return returnvalue;
		}

		double* _lb = new double[nV];
		returnvalue = Utilities::readFromFile( _lb, nV, lb_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _lb;
			return returnvalue;
		}

		double* _ub = new double[nV];
		returnvalue = Utilities::readFromFile( _ub, nV, ub_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _ub;
			return returnvalue;
		}

		double* _A = new double[nC*nV];
		returnvalue = Utilities::readFromFile( _A, nC*nV, A_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _A;
			return returnvalue;
		}
				
		double* _lbA = new double[nC];
		returnvalue = Utilities::readFromFile( _lbA, nC, lbA_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _lbA;
			return returnvalue;
		}

		double* _ubA = new double[nC];
		returnvalue = Utilities::readFromFile( _ubA, nC, ubA_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _ubA;
			return returnvalue;
		}

		double* _S1 = new double[nComp*nV];
		returnvalue = Utilities::readFromFile( _S1, nComp*nV, S1_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _S1;
			return returnvalue;
		}

		double* _S2 = new double[nComp*nV];
		returnvalue = Utilities::readFromFile( _S2, nComp*nV, S2_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _S2;
			return returnvalue;
		}

		double* _x0 = new double[nC + 2*nComp];
		returnvalue = Utilities::readFromFile( _x0, nC + 2*nComp, x0_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _S1;
			return returnvalue;
		}

		double* _y0 = new double[nC + 2*nComp];
		returnvalue = Utilities::readFromFile( _y0, nC + 2*nComp, y0_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _S2;
			return returnvalue;
		}

		// Fill vaues
		returnvalue = setH( _H );
		
		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;

		returnvalue = setG( _g );
		
		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
			
		returnvalue = setLB( _lb );
		
		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
			
		returnvalue = setUB( _ub );
		
		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
			
		returnvalue = setConstraints( _A, _S1, _S2, _lbA, _ubA );
		
		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;

		returnvalue = setInitialGuess( _x0, _y0 );
		
		if (returnvalue != SUCCESSFUL_RETURN)
			return returnvalue;
			
		return SUCCESSFUL_RETURN;
	}



	/*
	*	s e t C o n s t r a i n t s
	*/
	returnValue LCQProblem::setConstraints( 	const double* const A_new, const double* const S1_new, const double* const S2_new, 
												const double* const lbA_new, const double* const ubA_new )
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


	/*
	 *	 r u n S o l v e r
	 */
	returnValue LCQProblem::runSolver( ) {
		
		// Initialize variables
		initializeSolver();
		
		// Initialization strategy
		if (options.solveZeroPenaltyFirst) {
			memcpy(gk, g, nV*sizeof(double));

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

			// Print iteration
			printIteration( );

			// Terminate, update pen, or continue inner loop
			if (stationarityCheck()) {
				if (complementarityCheck()) {
					// Switch from penalized to LCQP duals
					transformDuals();

					// Determine C-,M-,S-Stationarity					
					algoStat = determineStationarityType();

					return SUCCESSFUL_RETURN;
				} else {
					updatePenalty();

					// Update iterate counters
					outerIter++;
					innerIter = 0;
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

	/*
	 *	 i n i t i a l i z e S o l v e r
	 */
	void LCQProblem::initializeSolver( ) {
		// Allocate vectors
		Qk = new double[nV*nV];
		gk = new double[nV];
		xnew = new double[nV];
		pk = new double[nV];
		statk = new double[nV];		
		
		// Initialize variables and counters
		alphak = 1;
		rho = options.initialComplementarityPenalty;		
		outerIter = 0;
		innerIter = 0;
		algoStat = algorithmStatus::PROBLEM_NOT_SOLVED;

		// Silent the subproblem solver
		qpOASES::Options qpOpts;
		if (options.printLvl < printLevel::VERBOSE)
			qpOpts.printLevel =  qpOASES::PrintLevel::PL_NONE;
		qp.setOptions( qpOpts );
	}


	/*
	 *	 s o l v e Q P S u b p r o b l e m
	 */
	returnValue LCQProblem::solveQPSubproblem(bool initialSolve) {
		// For now, only implementing for qpOASES solver
		qpOASES::returnValue ret;
		qpOASES::int_t nwsr = 10000;

		if (initialSolve) {
			ret = qp.init(H, g, A, lb, ub, lbA, ubA, nwsr, (double*)0, xk, yk);

			if (yk == 0)
				yk = new double[nV + nC + 2*nComp];

			hessianType = qp.getHessianType();
		} else {
			qp.setHessianType(hessianType);

			ret = qp.hotstart(gk, lb, ub, lbA, ubA, nwsr);
		}	

		if (ret != qpOASES::SUCCESSFUL_RETURN)
			return SUBPROBLEM_SOLVER_ERROR;

		// Update xnew, yk
		qp.getPrimalSolution(xnew);
		qp.getDualSolution(yk);

		// Update pk
		Utilities::WeightedVectorAdd(1, xnew, -1, xk, pk, nV);
		return SUCCESSFUL_RETURN;
	}

	
	/*
	 *	 s t a t i o n a r i t y C h e c k
	 */
	bool LCQProblem::stationarityCheck( ) {
		return Utilities::MaxAbs(statk, nV) < options.stationarityTolerance;
	}


	/*
	 *	 c o m p l e m e n t a r i t y C h e c k
	 */
	bool LCQProblem::complementarityCheck( ) {
		return Utilities::QuadraticFormProduct(C, xk, nV) < 2*options.complementarityTolerance;
	}


	/*
	 *	 t r a n s f o r m D u a l s
	 */
	void LCQProblem::transformDuals( ) {

		double* tmp = new double[nV];
		
		// y_S1 = y - rho*S2*xk
		Utilities::MatrixMultiplication(S2, xk, tmp, nV, nComp, 1);
		for (int i = 0; i < nComp; i++) {
			yk[nC + i] = yk[nC + i] - rho*tmp[i];
		}

		// y_S2 = y - rho*S1*xk
		Utilities::MatrixMultiplication(S1, xk, tmp, nComp, nV, 1);
		for (int i = 0; i < nComp; i++) {
			yk[nC + i] = yk[nC + i] - rho*tmp[i];
		}
	}


	/*
	 *	 d e t e r m i n e S t a t i o n a r i t y T y p e
	 */
	algorithmStatus LCQProblem::determineStationarityType( ) {
		//TODO: Need to get weakly active constraints first
		std::cout << "NOT YET IMPLEMENTED" << std::endl;
		
		return algorithmStatus::C_STATIONARY_SOLUTION;
	}


	/*
	 *	 u p d a t e P e n a l t y
	 */
	void LCQProblem::updatePenalty( ) {
		rho *= options.complementarityPenaltyUpdate;
	}


	/*
	 *	 g e t O p t i m a l S t e p L e n g t h
	 */
	void LCQProblem::getOptimalStepLength( ) {

		double qk = Utilities::QuadraticFormProduct(Qk, pk, nV);

		double* lk_tmp = new double[nV];
		Utilities::AffineLinearTransformation(1, Qk, xk, g, lk_tmp, nV, nV);

		double lk = Utilities::DotProduct(pk, lk_tmp, nV);
		
		alphak = 0;

		// Non convex case
		if (qk <= Utilities::EPS) {
			if (qk + lk < 0)
				alphak = 1;
		} else {
			// Descent + Convex
			if (lk < 0) {
				alphak = std::min(-lk/qk, 1.0);
			}
		}

		// 0-Step Length:
		if (Utilities::MaxAbs(pk, nV) < options.stationarityTolerance || complementarityCheck()) {
			alphak = 1;
		}
	}


	/*
	 *	 u p d a t e S t e p
	 */
	void LCQProblem::updateStep( ) {
		// Update penalty on rejected step
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
		}

		// Update Qk
		Utilities::WeightedMatrixAdd(1, H, rho, C, Qk, nV, nV);

		// stat = Qk*xk + g - A*yk
		// 1) Qk*xk
		Utilities::MatrixMultiplication(Qk, xk, statk, nV, nV, 1);

		// 2) -1*A*yk + g
		double* tmp_stat = new double[nV];
		Utilities::AffineLinearTransformation(-1, A, yk, g, tmp_stat, nV, nC + 2*nComp);
		Utilities::WeightedVectorAdd(1, tmp_stat, 1, statk, statk, nV);		
	}


	/*
	 *	 p e r t u r b G r a d i e n t
	 */
	void LCQProblem::perturbGradient( ) {
		
		int randNum;
		for (int i = 0; i < nV; i++) {
			// Random number -1, 0, 1
			randNum = (rand() % 3) - 1;

			gk[i] += randNum*Utilities::EPS;
		}
	}



	/*
	*	g e t P r i m a l S o l u t i o n
	*/
	algorithmStatus LCQProblem::getPrimalSolution( double* const xOpt ) const
	{
		if (algoStat != algorithmStatus::PROBLEM_NOT_SOLVED) {
			for (int i = 0; i < nV; i++)
				xOpt[i] = xk[i];
		}
		
		return algoStat;
	}


	/*
	*	g e t D u a l S o l u t i o n
	*/
	algorithmStatus LCQProblem::getDualSolution( double* const yOpt ) const
	{
		if (algoStat != algorithmStatus::PROBLEM_NOT_SOLVED) {
			for (int i = 0; i < nC + 2*nComp; i++)
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

		auto sep = " | ";

		// Print outer iterate
		printf("%6d", outerIter);

		// Print innter iterate
		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			printf("%s%*d", sep, 6, innerIter);

		printf("%s%10.3e%s%10.3e%s%10.3e", sep, Utilities::MaxAbs(statk, nV), sep, rho, sep, Utilities::MaxAbs(pk, nV));

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
			printf("%s%10.3e%s%*d", sep, alphak, sep, 6, qpIterk);
		}	

		printf(" \n");
	}

	/*
	 * 	 p r i n t H e a d e r 
	 */ 
	void LCQProblem::printHeader() 
	{
		printLine();

		auto sep = " | ";
		auto outer = " outer";
		auto inner = " inner";
		auto stat = "   stat   ";
		auto pen = "  penalty ";
		auto np = "  norm p  ";
		auto sl = " step len ";
		auto subIt = "sub it";

		printf("%s",outer);

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			printf("%s%s", sep, inner);

		printf("%s%s%s%s%s%s", sep, stat, sep, pen, sep, np);

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
		auto iSep = "------";
		auto dSep = "----------";
		auto node = "-+-";

		printf("%s", iSep);

		// Print innter iterate
		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			printf("%s%s", node, iSep);

		printf("%s%s%s%s%s%s", node, dSep, node, dSep, node, dSep);

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
			printf("%s%s%s%s", node, dSep, node, iSep);
		}	
		
		printf("-\n");
	}
}


	


/*
 *	end of file
 */