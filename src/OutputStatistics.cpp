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


#include "OutputStatistics.hpp"

namespace LCQPanther {

    OutputStatistics::OutputStatistics( ) { }


    OutputStatistics& OutputStatistics::operator=( const OutputStatistics& rhs )
    {
        iter_total = rhs.iter_total;
        iter_outer = rhs.iter_outer;
        subproblem_iter = rhs.subproblem_iter;
        rho_opt = rhs.rho_opt;
        status = rhs.status;

        return *this;
    }


    void OutputStatistics::reset( )
    {
        iter_total = 0;
        iter_outer = 0;
        subproblem_iter = 0;
        rho_opt = 0.0;
        status = PROBLEM_NOT_SOLVED;
    }


    ReturnValue OutputStatistics::updateIterTotal( int delta_iter )
    {
        if (delta_iter < 0) return INVALID_TOTAL_ITER_COUNT;

        iter_total += delta_iter;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateIterOuter( int delta_iter )
    {
        if (delta_iter < 0) return INVALID_TOTAL_OUTER_ITER;

        iter_outer += delta_iter;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateSubproblemIter( int delta_iter )
    {
        if (delta_iter < 0) return IVALID_SUBPROBLEM_ITER;

        subproblem_iter += delta_iter;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateRhoOpt( double _rho )
    {
        if (_rho <= 0) return INVALID_RHO_OPT;

        rho_opt = _rho;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue OutputStatistics::updateSolutionStatus( AlgorithmStatus _status )
    {
        status = _status;
        return SUCCESSFUL_RETURN;
    }


    int OutputStatistics::getIterTotal( ) const
    {
        return iter_total;
    }


    int OutputStatistics::getIterOuter( ) const
    {
        return iter_outer;
    }


    int OutputStatistics::getSubproblemIter( ) const
    {
        return subproblem_iter;
    }


    double OutputStatistics::getRhoOpt( ) const
    {
        return rho_opt;
    }


    AlgorithmStatus OutputStatistics::getSolutionStatus( ) const
    {
        return status;
    }
}
