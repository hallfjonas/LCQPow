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


namespace lcqpOASES {
	/*****************************************************************************
	 *  P U B L I C                                                              *
	 *****************************************************************************/
	/*
	*	g e t N C o m p
	*/
	inline int LCQProblem::getNV( ) const
	{
		return nV;
	}


	/*
	*	g e t N C o m p
	*/
	inline int LCQProblem::getNC( ) const
	{
		return nC;
	}


	/*
	*	g e t N C o m p
	*/
	inline int LCQProblem::getNComp( ) const
	{
		return nComp;
	}


	/*
	*	s e t H
	*/
	inline returnValue LCQProblem::setH( const double* const H_new )
	{
		uint nV = (uint)getNV( );

		if (nV <= 0)
			return returnValue::LCQPOBJECT_NOT_SETUP;

		H = new double[nV*nV];
		memcpy( H, H_new, nV*nV*sizeof(double) );

		return returnValue::SUCCESSFUL_RETURN;
	}


	/*
	*	s e t G
	*/
	inline returnValue LCQProblem::setG( const double* const g_new )
	{
		uint nV = (uint)getNV( );

		if ( nV == 0 )
			return returnValue::LCQPOBJECT_NOT_SETUP;

		if ( g_new == 0 )
			return returnValue::ILLEGAL_ARGUMENT;

		g = new double[nV];
		memcpy( g, g_new, nV*sizeof(double) );

		return returnValue::SUCCESSFUL_RETURN;
	}


	/*
	*	s e t L B
	*/
	inline returnValue LCQProblem::setLB( const double* const lb_new )
	{
		uint i;
		uint nV = (uint)getNV( );

		if ( nV == 0 )
			return returnValue::LCQPOBJECT_NOT_SETUP;

		lb = new double[nV];

		if ( lb_new != 0 )
		{
			memcpy( lb,lb_new,nV*sizeof(double) );
		}
		else
		{
			/* if no lower bounds are specified, set them to -infinity */
			for( i = 0; i < nV; i++ )
				lb[i] = -INFINITY;
		}

		return returnValue::SUCCESSFUL_RETURN;
	}


	/*
	*	s e t L B
	*/
	inline returnValue LCQProblem::setLB( int number, double value )
	{
		int nV = getNV( );

		if ( nV == 0 )
			return returnValue::LCQPOBJECT_NOT_SETUP;

		if ( ( number >= 0 ) && ( number < nV ) )
		{
			lb[number] = value;
			return returnValue::SUCCESSFUL_RETURN;
		}
		else
		{
			return returnValue::INDEX_OUT_OF_BOUNDS;
		}
	}


	/*
	*	s e t U B
	*/
	inline returnValue LCQProblem::setUB( const double* const ub_new )
	{
		uint i;
		uint nV = (uint)getNV( );

		if ( nV == 0 )
			return returnValue::LCQPOBJECT_NOT_SETUP;

		ub = new double[nV];

		if ( ub_new != 0 )
		{
			memcpy( ub,ub_new,nV*sizeof(double) );
		}
		else
		{
			/* if no upper bounds are specified, set them to infinity */
			for( i=0; i<nV; ++i )
				ub[i] = INFINITY;
		}

		return returnValue::SUCCESSFUL_RETURN;
	}


	/*
	*	s e t U B
	*/
	inline returnValue LCQProblem::setUB( int number, double value )
	{
		int nV = getNV( );

		if ( nV == 0 )
			return returnValue::LCQPOBJECT_NOT_SETUP;

		if ( ( number >= 0 ) && ( number < nV ) )
		{
			ub[number] = value;

			return returnValue::SUCCESSFUL_RETURN;
		}
		else
		{
			return returnValue::INDEX_OUT_OF_BOUNDS;
		}
	}

	/*
	 *	s e t O p t i o n s
	 */
	inline void LCQProblem::setOptions( const Options& _options )
	{
		options = _options;
		options.ensureConsistency( );
	}
}

/*
 *	end of file
 */
