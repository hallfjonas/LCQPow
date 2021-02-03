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
#include <iomanip>
#include <qpOASES.hpp>

using qpOASES::QProblem;
using std::cout;
using std::setw;
using std::string;

namespace lcqpOASES {

	/*****************************************************************************
	 *  P U B L I C                                                              *
	 *****************************************************************************/

	/*
	*	Q P r o b l e m
	*/
	LCQProblem::LCQProblem( int _nV, int _nC, int _nComp )
	{
		int i;

		/* consistency checks */
		if ( _nV <= 0 || _nComp <= 0)
		{
			throw( ILLEGAL_ARGUMENT );
		}

		if ( _nC < 0 )
		{
			_nC = 0;
			throw( ILLEGAL_ARGUMENT );
		}

		nV = _nV;
		nC = _nC;
		nComp = _nComp;

		QProblem lp(nV, nC);
		
		int nc = lp.getNAC();

		qp = lp;
	}


	/*
	*	s o l v e
	*/
	returnValue LCQProblem::solve(	const double* const _H, const double* const _g, const double* const _A,
									const double* const _lb, const double* const _ub,
									const double* const _lbA, const double* const _ubA,
									const double* const _S1, const double* const _S2,
									double* const cputime, const double* const xOpt, const double* const yOpt
									)
	{
		setupLCQPdata(_H, _g, _A, _lb, _ub, _lbA, _ubA, _S1, _S2);

		int nV = getNV();
		int nC = getNC();
		int nComp = getNComp();

		qpOASES::int_t nwsr = 1000;

		for (int i = 0; i < 20; i++) {
			for (int k = 0; k < 30; k++) {
				printIteration(i, 10, 100, 0.3, 0.00003, k, 2, 60);		
			}
		}
		

		if (options.solveZeroPenaltyFirst)
			qpOASES::returnValue ret = qp.init(H, g, A, lb, ub, lbA, ubA, nwsr);		

		// TODO: Write solver
		return returnValue::NOT_YET_IMPLEMENTED;	
	}


	/*
	*	i n i t
	*/
	returnValue LCQProblem::solve(	const char* const H_file, const char* const g_file, const char* const A_file,
									const char* const lb_file, const char* const ub_file,
									const char* const lbA_file, const char* const ubA_file,
									const char* const S1_file, const char* const S2_file,
									int& nWSR, double* const cputime,
									const double* const xOpt, const double* const yOpt
									)
	{
		// Internally set up data
		setupLCQPdata(H_file, g_file, A_file, lb_file, ub_file, lbA_file, ubA_file, S1_file, S2_file);

		// TODO: Write this
		return returnValue::SUCCESSFUL_RETURN;	
	}


	/*
	*	s e t C o n s t r a i n t s
	*/
	returnValue LCQProblem::setConstraints( 	const double* const A_new, const double* const S1_new, const double* const S2_new, 
												const double* const lbA_new, const double* const ubA_new )
	{
		int j;
		int nV = getNV( );
		int nC = getNC( );
		int nComp = getNComp( );
		
		if ( nV == 0 || nComp == 0 )
			return returnValue::LCQPOBJECT_NOT_SETUP;

		if ( A_new == 0 && nC > 0)
			return returnValue::ILLEGAL_ARGUMENT;

		// Set up new constraint matrix (A; L; R)
		A = new double[(nC + 2*nV)*nV];

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
			return returnValue::ILLEGAL_ARGUMENT;

		S1 = new double[nComp*nV];
		S2 = new double[nComp*nV];

		for (int i = 0; i < nComp*nV; i++) {
			S1[i] = S1_new[i];
			S2[i] = S2_new[i];
		}

		Utilities utils;
		double* C = new double[nV*nV];
		utils.MatrixSymmetrizationProduct(S1, S2, C, nComp, nV);
		
		return returnValue::SUCCESSFUL_RETURN;
	}


	/*
	*	g e t O b j V a l
	*/
	double LCQProblem::getObjVal( ) const
	{
		return qp.getObjVal( );
	}


	/*
	*	g e t O b j V a l
	*/
	double LCQProblem::getObjVal( const double* const _x ) const
	{
		return qp.getObjVal( _x );
	}


	/*
	*	g e t P r i m a l S o l u t i o n
	*/
	returnValue LCQProblem::getPrimalSolution( double* const xOpt ) const
	{
		qpOASES::returnValue qpRet = qp.getPrimalSolution( xOpt );

		if (qpRet == qpOASES::returnValue::SUCCESSFUL_RETURN)
			return returnValue::SUCCESSFUL_RETURN;

		return returnValue::SUBPROBLEM_SOLVER_ERROR;
	}


	/*
	*	g e t D u a l S o l u t i o n
	*/
	returnValue LCQProblem::getDualSolution( double* const yOpt ) const
	{
		// TODO: Implement transformation
		return getDualSolution( yOpt );
	}


	/*
	 *	 p r i n t I t e r a t i o n
	 */
	void LCQProblem::printIteration(	int outerIter,							/**< Number of current outer iteration. */
										double stationarityValue,				/**< Stationarity value of outer loop problem. */
										double complementarityValue,			/**< Current complementarity value. */
										double penaltyValue,					/**< Current penalty value. */
										double normStep,						/**< Euclidean distance to last iterate. */
										int innerIter,							/**< Number of current inner iteration. */
										double optimalStepLength,   			/**< Step length of current iteration computed via optimal step length approach. */
										int qpIter		 						/**< Number of iterations performed by subproblem solver. */
										) 
	{
		if (options.printLvl == printLevel::NONE)
			return;

		// Print header every 10 iters
		bool headerInner = (options.printLvl >= printLevel::INNER_LOOP_ITERATES && innerIter % 10 == 0);
		bool headerOuter = (options.printLvl == printLevel::OUTER_LOOP_ITERATES && outerIter % 10 == 0);
		if (headerInner || headerOuter)
			printHeader();	

		int iL = printIntLength;
		int dL = printDoubleLength;

		// Print outer iterate
		cout << setw(iL) << outerIter;

		// Print innter iterate
		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			cout << " | " << setw(iL) << innerIter;

		cout << " | " << setw(dL) << stationarityValue;

		cout << " | " << setw(dL) << penaltyValue;

		cout << " | " << setw(dL) << normStep;

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
			cout << " | " << setw(dL) << optimalStepLength;
			cout << " | " << setw(iL) << qpIter;
		}	

		cout << std::endl;		
	}


	void LCQProblem::printHeader() 
	{
		printLine();

		int iL = printIntLength;
		int dL = printDoubleLength;

		cout << setw(iL) << "outer";

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			cout << " | " << setw(iL) << "inner";

		cout << " | " << setw(dL) << "stat";
		cout << " | " << setw(dL) << "penalty";
		cout << " | " << setw(dL) << "norm p";

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
			cout << " | " << setw(dL) << "step len";
			cout << " | " << setw(iL) << "sub it";
		}	
		cout << std::endl;		

		printLine();
	}

	void LCQProblem::printLine() 
	{
	
		int iL = printIntLength;
		int dL = printDoubleLength;

		// Print line
		cout << setw(iL) << string(iL, '-');

		// Print innter iterate
		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES)
			cout << "-+-" << setw(iL) << string(iL, '-');

		cout << "-+-" << setw(dL) << string(dL, '-');
		cout << "-+-" << setw(dL) << string(dL, '-');
		cout << "-+-" << setw(dL) << string(dL, '-');

		if (options.printLvl >= printLevel::INNER_LOOP_ITERATES) {
			cout << "-+-" << setw(dL) << string(dL, '-');
			cout << "-+-" << setw(iL) << string(iL, '-');
		}	
		cout << std::endl;		
	}

	/*
	*	s e t u p Q P d a t a
	*/
	returnValue LCQProblem::setupLCQPdata(	const double* const _H, const double* const _g, const double* const _A,
											const double* const _lb, const double* const _ub,
											const double* const _lbA, const double* const _ubA,
											const double* const _S1, const double* const _S2
											)
	{
		setH( _H );
		setG( _g );
		setLB( _lb );
		setUB( _ub );
		setConstraints( _A, _S1, S2, _lbA, _ubA );	

		return returnValue::SUCCESSFUL_RETURN;
	}


	/*
	*	s e t u p Q P d a t a F r o m F i l e
	*/
	returnValue LCQProblem::setupLCQPdata(	const char* const H_file, const char* const g_file, const char* const A_file,
											const char* const lb_file, const char* const ub_file,
											const char* const lbA_file, const char* const ubA_file,
											const char* const S1_file, const char* const S2_file
											)
	{
		int i;
		int nV = getNV( );
		int nC = getNC( );

		returnValue returnvalue;
		Utilities utils;

		double* _H = new double[nV*nV];
		returnvalue = utils.readFromFile( _H, nV*nV, H_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _H;
			throw( returnvalue );
		}
			

		double* _g = new double[nV];
		returnvalue = utils.readFromFile( _g, nV, g_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _g;
			throw( returnvalue );
		}

		double* _lb = new double[nV];
		returnvalue = utils.readFromFile( _lb, nV, lb_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _lb;
			throw( returnvalue );
		}

		double* _ub = new double[nV];
		returnvalue = utils.readFromFile( _ub, nV, ub_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _ub;
			throw( returnvalue );
		}

		double* _A = new double[nC*nV];
		returnvalue = utils.readFromFile( _A, nC*nV, A_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _A;
			throw( returnvalue );
		}
				
		double* _lbA = new double[nC];
		returnvalue = utils.readFromFile( _lbA, nC, lbA_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _lbA;
			throw( returnvalue );
		}

		double* _ubA = new double[nC];
		returnvalue = utils.readFromFile( _ubA, nC, ubA_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _ubA;
			throw( returnvalue );
		}

		double* _S1 = new double[nComp*nV];
		returnvalue = utils.readFromFile( _S1, nComp*nV, S1_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _S1;
			throw( returnvalue );
		}

		double* _S2 = new double[nComp*nV];
		returnvalue = utils.readFromFile( _S2, nComp*nV, S2_file );
		if ( returnvalue != SUCCESSFUL_RETURN ) {
			delete _S2;
			throw( returnvalue );
		}


		// Fill vaues
		setH( _H );
		setG( _g );
		setLB( _lb );
		setUB( _ub );
		setConstraints( _A, _S1, _S2, _lbA, _ubA );

		return returnValue::SUCCESSFUL_RETURN;
	}
}


/*
 *	end of file
 */