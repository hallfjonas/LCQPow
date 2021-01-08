/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file src/LCQProblem.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the LCQProblem class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */


#include <qpOASES/LCQProblem.hpp>
#include <qpOASES/ConstraintProduct.hpp>
#include <qpOASES/LapackBlasReplacement.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

/*
 *	Q P r o b l e m
 */
LCQProblem::LCQProblem( int_t _nV, int_t _nC, int_t _nComp, HessianType _hessianType, BooleanType allocDenseMats )
{
	int_t i;

	/* consistency checks */
	if ( _nV <= 0 || _nComp <= 0)
	{
		_nV = 1; _nComp = 1;
		THROWERROR( RET_INVALID_ARGUMENTS );
	}

	if ( _nC < 0 )
	{
		_nC = 0;
		THROWERROR( RET_INVALID_ARGUMENTS );
	}

	nV = _nV;
	nC = _nC;
	nComp = _nComp;

	QProblem _qp(_nV, _nC, _hessianType, allocDenseMats);
	qp = _qp;
}


/*
 *	i n i t
 */
returnValue LCQProblem::solve(	SymmetricMatrix *_H, const real_t* const _g, Matrix *_A,
								const real_t* const _lb, const real_t* const _ub,
								const real_t* const _lbA, const real_t* const _ubA,
								Matrix *_S1, Matrix *_S2, 
								int_t& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
								const real_t* const _R
								)
{
	setupLCQPdata(_H, _g, _A, _lb, _ub, _lbA, _ubA, _S1, _S2);

	int_t nV = getNV();
	int_t nC = getNC();
	int_t nComp = getNComp();

	// TODO 
}


/*
 *	i n i t
 */
returnValue LCQProblem::solve(	const real_t* const _H, const real_t* const _g, const real_t* const _A,
								const real_t* const _lb, const real_t* const _ub,
								const real_t* const _lbA, const real_t* const _ubA,
								const real_t* const _S1, const real_t* const _S2,
								int_t& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
								const real_t* const _R
								)
{
	setupLCQPdata(_H, _g, _A, _lb, _ub, _lbA, _ubA, _S1, _S2);

	int_t nV = getNV();
	int_t nC = getNC();
	int_t nComp = getNComp();

	// TODO: Write solver
	return RET_NOT_YET_IMPLEMENTED;	
}


/*
 *	i n i t
 */
returnValue LCQProblem::solve(	const char* const H_file, const char* const g_file, const char* const A_file,
								const char* const lb_file, const char* const ub_file,
								const char* const lbA_file, const char* const ubA_file,
								const char* const S1_file, const char* const S2_file,
								int_t& nWSR, real_t* const cputime,
								const real_t* const xOpt, const real_t* const yOpt,
								const Bounds* const guessedBounds, const Constraints* const guessedConstraints,
								const char* const R_file
								)
{
	// Internally set up data
	setupLCQPdata(H_file, g_file, A_file, lb_file, ub_file, lbA_file, ubA_file, S1_file, S2_file);

	// TODO: Write this
	return SUCCESSFUL_RETURN;	
}


/*
 *	s e t P r i n t L e v e l
 */
returnValue LCQProblem::setPrintLevel( PrintLevel _printLevel )
{
	#ifndef __SUPPRESSANYOUTPUT__
		#ifndef __MATLAB__
			if ( ( options.printLevel == PL_HIGH ) && ( options.printLevel != _printLevel ) )
				THROWINFO( RET_PRINTLEVEL_CHANGED );
		#endif /* __MATLAB__ */
		options.printLevel = _printLevel;
	#else
	options.printLevel = PL_NONE;
	#endif /* __SUPPRESSANYOUTPUT__ */

	/* update message handler preferences */
 	switch ( options.printLevel )
 	{
 		case PL_NONE:
 			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		case PL_TABULAR:
		case PL_LOW:
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_HIDDEN );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		case PL_DEBUG_ITER:
		case PL_MEDIUM:
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_HIDDEN );
			break;

		default: /* PL_HIGH */
			getGlobalMessageHandler( )->setErrorVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setWarningVisibilityStatus( VS_VISIBLE );
			getGlobalMessageHandler( )->setInfoVisibilityStatus( VS_VISIBLE );
			break;
 	}

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t O b j V a l
 */
real_t LCQProblem::getObjVal( ) const
{
	return qp.getObjVal( );
}


/*
 *	g e t O b j V a l
 */
real_t LCQProblem::getObjVal( const real_t* const _x ) const
{
	return qp.getObjVal( _x );
}


/*
 *	g e t P r i m a l S o l u t i o n
 */
returnValue LCQProblem::getPrimalSolution( real_t* const xOpt ) const
{
	return qp.getPrimalSolution( xOpt );
}


/*
 *	g e t D u a l S o l u t i o n
 */
returnValue LCQProblem::getDualSolution( real_t* const yOpt ) const
{
	// TODO: Implement transformation
	return getDualSolution( yOpt );
}


/*
 *	s e t u p Q P d a t a
 */
returnValue LCQProblem::setupLCQPdata(	SymmetricMatrix *_H, const real_t* const _g, Matrix *_A,
										const real_t* const _lb, const real_t* const _ub,
										const real_t* const _lbA, const real_t* const _ubA,
										Matrix *_S1, Matrix *_S2
										)
{
	setH( _H );
	setG( _g );
	setLB( _lb );
	setUB( _ub );
	setComplementarities(_S1, _S2);
	setConstraints(_A, _lbA, _ubA);

	return SUCCESSFUL_RETURN;
}

/*
 *	s e t u p Q P d a t a
 */
returnValue LCQProblem::setupLCQPdata(	const real_t* const _H, const real_t* const _g, const real_t* const _A,
										const real_t* const _lb, const real_t* const _ub,
										const real_t* const _lbA, const real_t* const _ubA,
										const real_t* const _S1, const real_t* const _S2
										)
{
	setH( _H );
	setG( _g );
	setLB( _lb );
	setUB( _ub );
	setComplementarities(_S1, _S2);
	setConstraints( _A, _lbA, _ubA );	

	return SUCCESSFUL_RETURN;
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
	int_t i;
	int_t nV = getNV( );
	int_t nC = getNC( );

	returnValue returnvalue;

	real_t* _H = new real_t[nV * nV];
	returnvalue = readFromFile( _H, nV, nV, H_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		delete[] _H;
		return THROWERROR( returnvalue );
	}

	returnvalue = readFromFile( g, nV, g_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );

	returnvalue = readFromFile( lb, nV, lb_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );

	returnvalue = readFromFile( ub, nV, ub_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );

	real_t* _S1 = new real_t[nComp * nV];
	returnvalue = readFromFile( _S1, nComp, nV, S1_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		delete[] _S1;
		return THROWERROR( returnvalue );
	}

	real_t* _S2 = new real_t[nComp * nV];
	returnvalue = readFromFile( _S2, nComp, nV, S2_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		delete[] _S2;
		return THROWERROR( returnvalue );
	}

	real_t* _lbA = new real_t[nC];
	returnvalue = readFromFile( _lbA, nC, lbA_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		delete[] _lbA;
		return THROWERROR( returnvalue );
	}

	real_t* _ubA = new real_t[nC];
	returnvalue = readFromFile( _ubA, nC, ubA_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		delete[] _ubA;
		return THROWERROR( returnvalue );
	}

	real_t* _A = new real_t[nC * nV];
	returnvalue = readFromFile( _A, nC, nV, A_file );
	if ( returnvalue != SUCCESSFUL_RETURN )
	{
		delete[] _A;
		return THROWERROR( returnvalue );
	}
	
	// Fill vaues
	setH( _H );
	setComplementarities(_S1, _S2);
	setConstraints( _A, _lbA, _ubA );
		
	H->doFreeMemory( );

	return SUCCESSFUL_RETURN;
}



returnValue LCQProblem::setC()
{
	int_t nV = getNV();
	int_t nComp = getNComp();

	real_t* C_new = new real_t[nV*nV];
	
	real_t* tmpS1 = S1->full();
	real_t* tmpS2 = S2->full();

	// tmp = S1'*S2
	S1->transTimes(nV, 1.0, tmpS2, nComp, 0.0, C_new, nV);

	// tmp = S2'*S1 + tmp 
	S2->transTimes(nV, 1.0, tmpS1, nComp, 1.0, C_new, nV);

	SymDenseMat* dC;
	C = dC = new SymDenseMat( nV, nV, nV, (real_t*) C_new );
}


returnValue LCQProblem::setConstraints(	Matrix* A_new, const real_t* const lbA_new, const real_t* const ubA_new)
{
	return RET_NOT_YET_IMPLEMENTED;	
}

returnValue LCQProblem::setConstraints(	const real_t* const A_new, const real_t* const lbA_new, const real_t* const ubA_new)
{

	int_t nV = getNV( );
	int_t nC = getNC( );
	int_t nComp = getNComp( );

	real_t* S1_full = S1->full();
	real_t* S2_full = S2->full();

	real_t* _A = new real_t[(nC + 2*nComp) * nV];
	lbA = new real_t[nC + 2*nComp];
	ubA = new real_t[nC + 2*nComp];

	// First, fill constraint matrix
	for (int i = 0; i < nC; i++) {
		lbA[i] = lbA_new[i];
		ubA[i] = ubA_new[i];

		for (int j = 0; j < nV; j++)
			_A[i*nV + j] = A_new[i*nV + j];	
	}

	// Then, fill values according to complementarity matrices
	for (int i = 0; i < nComp; i++) {
		int_t rowS1 = nC + i;
		int_t rowS2 = nC + nComp + i;
	
		// Fill lower bounds
		lbA[rowS1] = 0;
		lbA[rowS2] = 0;

		// Fill upper bounds
		ubA[rowS1] = INFINITY;
		ubA[rowS2] = INFINITY;

		// Fill matrix
		for (int j = 0; j < nV; j++) {
			_A[rowS1*nV + j] = S1_full[i*nV + j];	
			_A[rowS2*nV + j] = S2_full[i*nV + j];
		}				
	}

	setA( _A, nC + 2*nComp );

	return SUCCESSFUL_RETURN;	
}

END_NAMESPACE_QPOASES


/*
 *	end of file
 */