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
 *	\file include/qpOASES/LCQProblem.ipp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of inlined member functions of the LCQProblem class which
 *	is able to use the newly developed online active set strategy for
 *	parametric quadratic programming.
 */

BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/
/*
 *	g e t N C o m p
 */
inline int_t LCQProblem::getNV( ) const
{
	return nV;
}

/*
 *	g e t N C o m p
 */
inline int_t LCQProblem::getNC( ) const
{
	return nC;
}

/*
 *	g e t N C o m p
 */
inline int_t LCQProblem::getNComp( ) const
{
	return nComp;
}

/*
 *	s e t H
 */
inline returnValue LCQProblem::setH( SymmetricMatrix* H_new )
{
	if ( ( freeHessian == BT_TRUE ) && ( H != 0 ) )
	{
		delete H;
		H = 0;
	}

	H = H_new;
	freeHessian = BT_FALSE;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t H
 */
inline returnValue LCQProblem::setH( const real_t* const H_new )
{
	int_t nV = getNV();
	SymDenseMat* dH;

	/* Not allowing 0 hessians here */
	if ( H_new == 0 )
	{
		return THROWERROR( RET_INVALID_ARGUMENTS );
	}
	
	if ( ( freeHessian == BT_TRUE ) && ( H != 0 ) )
		delete H;

	H = dH = new SymDenseMat( nV, nV, nV, (real_t*) H_new );
	freeHessian = BT_TRUE;


	return SUCCESSFUL_RETURN;
}

/*
 *	s e t G
 */
inline returnValue LCQProblem::setG( const real_t* const g_new )
{
	uint_t nV = (uint_t)getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( g_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	g = new real_t[nV];
	memcpy( g, g_new, nV*sizeof(real_t) );

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t L B
 */
inline returnValue LCQProblem::setLB( const real_t* const lb_new )
{
	uint_t i;
	uint_t nV = (uint_t)getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	lb = new real_t[nV];

	if ( lb_new != 0 )
	{
		memcpy( lb,lb_new,nV*sizeof(real_t) );
	}
	else
	{
		/* if no lower bounds are specified, set them to -infinity */
		for( i=0; i<nV; ++i )
			lb[i] = -INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t L B
 */
inline returnValue LCQProblem::setLB( int_t number, real_t value )
{
	int_t nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ( number >= 0 ) && ( number < nV ) )
	{
		lb[number] = value;
		return SUCCESSFUL_RETURN;
	}
	else
	{
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );
	}
}


/*
 *	s e t U B
 */
inline returnValue LCQProblem::setUB( const real_t* const ub_new )
{
	uint_t i;
	uint_t nV = (uint_t)getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	ub = new real_t[nV];

	if ( ub_new != 0 )
	{
		memcpy( ub,ub_new,nV*sizeof(real_t) );
	}
	else
	{
		/* if no upper bounds are specified, set them to infinity */
		for( i=0; i<nV; ++i )
			ub[i] = INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t U B
 */
inline returnValue LCQProblem::setUB( int_t number, real_t value )
{
	int_t nV = getNV( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ( number >= 0 ) && ( number < nV ) )
	{
		ub[number] = value;

		return SUCCESSFUL_RETURN;
	}
	else
	{
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );
	}
}


/*
 *	s e t A
 */
inline returnValue LCQProblem::setA( Matrix *A_new )
{
	int_t j;
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( A_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* Set constraint matrix AND update member AX. */
	if ( ( freeConstraintMatrix == BT_TRUE ) && ( A != 0 ) )
	{
		delete A;
		A = 0;
	}
	A = A_new;
	freeConstraintMatrix = BT_FALSE;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t A
 */
inline returnValue LCQProblem::setA( const real_t* const A_new, const int_t nRows )
{
	int_t j;
	int_t nV = getNV( );
	DenseMatrix* dA;

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( A_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	/* Set constraint matrix AND update member AX. */
	if ( ( freeConstraintMatrix == BT_TRUE ) && ( A != 0 ) )
	{
		delete A;
		A = 0;
	}

	A = dA = new DenseMatrix(nRows, nV, nV, (real_t*) A_new);
	freeConstraintMatrix = BT_TRUE;

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t L B A
 */
inline returnValue LCQProblem::setLBA( const real_t* const lbA_new )
{
	uint_t i;
	uint_t nV = (uint_t)getNV( );
	uint_t nC = (uint_t)getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	lbA = new real_t[nC];

	if ( lbA_new != 0 )
	{
		memcpy( lbA,lbA_new,nC*sizeof(real_t) );
	}
	else
	{
		/* if no lower constraints' bounds are specified, set them to -infinity */
		for( i=0; i<nC; ++i )
			lbA[i] = -INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t L B A
 */
inline returnValue LCQProblem::setLBA( int_t number, real_t value )
{
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ( number >= 0 ) && ( number < nC ) )
	{
		lbA[number] = value;
		return SUCCESSFUL_RETURN;
	}
	else
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );
}


/*
 *	s e t U B A
 */
inline returnValue LCQProblem::setUBA( const real_t* const ubA_new )
{
	uint_t i;
	uint_t nV = (uint_t)getNV( );
	uint_t nC = (uint_t)getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	ubA = new real_t[nC];

	if ( ubA_new != 0 )
	{
		memcpy( ubA,ubA_new,nC*sizeof(real_t) );
	}
	else
	{
		/* if no upper constraints' bounds are specified, set them to infinity */
		for( i=0; i<nC; ++i )
			ubA[i] = INFTY;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t U B A
 */
inline returnValue LCQProblem::setUBA( int_t number, real_t value )
{
	int_t nV = getNV( );
	int_t nC = getNC( );

	if ( nV == 0 )
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( ( number >= 0 ) && ( number < nC ) )
	{
		ubA[number] = value;
		return SUCCESSFUL_RETURN;
	}
	else
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );
}


/*
 *	s e t C o m p l e m e n t a r i t i e s
 */
inline returnValue LCQProblem::setComplementarities( Matrix* S1_new, Matrix* S2_new )
{
	// Free space if required
	if ( freeComplementarityMatrix == BT_TRUE )
	{
		if ( S1 != 0 ) {
			delete S1;
			S1 = 0;
		}

		if ( S2 != 0 ) {
			delete S2;
			S2 = 0;
		}
	}

	// Assign matrices
	S1 = S1_new;
	S2 = S2_new;

	// Set complementarity matrix
	setC();

	return SUCCESSFUL_RETURN;
}

/*
 *	s e t C o m p l e m e n t a r i t i e s
 */
inline returnValue LCQProblem::setComplementarities( const real_t* const S1_new,  const real_t* const S2_new )
{
	int_t j;
	int_t nV = getNV( );
	int_t nComp = getNComp( );
	DenseMatrix* dS1;
	DenseMatrix* dS2;

	if ( nV == 0 || nComp == 0)
		return THROWERROR( RET_QPOBJECT_NOT_SETUP );

	if ( S1_new == 0 || S2_new == 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	// Free space if required
	if ( freeComplementarityMatrix == BT_TRUE )
	{
		if ( S1 != 0 ) {
			delete S1;
			S1 = 0;
		}

		if ( S2 != 0 ) {
			delete S2;
			S2 = 0;
		}
	}

	// Assign values
	S1 = dS1 = new DenseMatrix(nComp, nV, nV, (real_t*) S1_new);
	S2 = dS2 = new DenseMatrix(nComp, nV, nV, (real_t*) S2_new);
	freeComplementarityMatrix = BT_TRUE;

	// Set complementarity matrix
	setC();

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t O p t i o n s
 */
inline returnValue LCQProblem::setOptions(	const Options& _options
											)
{
	options = _options;
	options.ensureConsistency( );

	setPrintLevel( options.printLevel );

	return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
