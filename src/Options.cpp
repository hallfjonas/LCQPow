/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2022 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "Options.hpp"
#include "MessageHandler.hpp"

namespace LCQPow {

    Options::Options( )
    {
        setToDefault( );
    }


    Options::Options( const Options& rhs)
    {
        copy( rhs );
    }


    Options::~Options( )
    { }


    Options& Options::operator=( const Options& rhs )
    {
        if ( this != &rhs )
        {
            copy( rhs );
        }

        return *this;
    }


    void Options::copy( const Options& rhs ) {
        stationarityTolerance = rhs.stationarityTolerance;
        complementarityTolerance = rhs.complementarityTolerance;
        initialPenaltyParameter = rhs.initialPenaltyParameter;
        penaltyUpdateFactor = rhs.penaltyUpdateFactor;
        solveZeroPenaltyFirst = rhs.solveZeroPenaltyFirst;
        maxIterations = rhs.maxIterations;
        maxRho = rhs.maxRho;
        nDynamicPenalty = rhs.nDynamicPenalty;
        etaDynamicPenalty = rhs.etaDynamicPenalty;
        printLevel = rhs.printLevel;
        storeSteps = rhs.storeSteps;
        qpSolver = rhs.qpSolver;
    }


    double Options::getStationarityTolerance( ) {
        return stationarityTolerance;
    }


    ReturnValue Options::setStationarityTolerance( double val ) {
        if (val <= Utilities::EPS)
            return (MessageHandler::PrintMessage(INVALID_STATIONARITY_TOLERANCE,WARNING) );

        stationarityTolerance = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    double Options::getComplementarityTolerance( ) {
        return complementarityTolerance;
    }


    ReturnValue Options::setComplementarityTolerance( double val ) {
        if (val <= Utilities::EPS)
            return (MessageHandler::PrintMessage(INVALID_COMPLEMENTARITY_TOLERANCE,WARNING) ) ;

        complementarityTolerance = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    double Options::getInitialPenaltyParameter( ) {
        return initialPenaltyParameter;
    }


    ReturnValue Options::setInitialPenaltyParameter( double val ) {
        if (val <= Utilities::ZERO)
            return (MessageHandler::PrintMessage(INVALID_INITIAL_PENALTY_VALUE,WARNING) );

        initialPenaltyParameter = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    double Options::getPenaltyUpdateFactor( ) {
        return penaltyUpdateFactor;
    }


    ReturnValue Options::setPenaltyUpdateFactor( double val ) {
        if (val <= 1)
            return (MessageHandler::PrintMessage(INVALID_PENALTY_UPDATE_VALUE,WARNING) ) ;

        penaltyUpdateFactor = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    bool Options::getSolveZeroPenaltyFirst( ) {
        return solveZeroPenaltyFirst;
    }


    ReturnValue Options::setSolveZeroPenaltyFirst( bool val ) {
        solveZeroPenaltyFirst = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    int Options::getMaxIterations( ) {
        return maxIterations;
    }


    ReturnValue Options::setMaxIterations( int val ) {
        if (val <= 0)
            return (MessageHandler::PrintMessage(INVALID_MAX_ITERATIONS_VALUE,WARNING) );

        maxIterations = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    double Options::getMaxRho( ) {
        return maxRho;
    }


    ReturnValue Options::setMaxRho( double val ) {
        if (val <= 0)
            return (MessageHandler::PrintMessage(INVALID_MAX_RHO_VALUE,WARNING) );

        maxRho = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }

    int Options::getNDynamicPenalty( ) {
        return nDynamicPenalty;
    }


    ReturnValue Options::setNDynamicPenalty( int val ) {
        nDynamicPenalty = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    double Options::getEtaDynamicPenalty( ) {
        return etaDynamicPenalty;
    }


    ReturnValue Options::setEtaDynamicPenalty( double val ) {
        if (val <= Utilities::EPS || val >= 1)
            return (MessageHandler::PrintMessage(INVALID_ETA_VALUE,WARNING) );

        etaDynamicPenalty = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    PrintLevel Options::getPrintLevel( ) {
        return printLevel;
    }


    ReturnValue Options::setPrintLevel( PrintLevel val ) {
        printLevel = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    ReturnValue Options::setPrintLevel( int val ) {

        if (val < PrintLevel::NONE || val > PrintLevel::SUBPROBLEM_SOLVER_ITERATES)
            return (MessageHandler::PrintMessage(INVALID_PRINT_LEVEL_VALUE,WARNING) );

        printLevel = (PrintLevel)val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    bool Options::getStoreSteps( ) {
        return storeSteps;
    }


    ReturnValue Options::setStoreSteps( bool val ) {
        storeSteps = val;
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    QPSolver Options::getQPSolver( ) {
        return qpSolver;
    }


    ReturnValue Options::setQPSolver( QPSolver val ) {
        qpSolver = val;
        return SUCCESSFUL_RETURN;
    }


    ReturnValue Options::setQPSolver( int val ) {
        if (val < QPSolver::QPOASES_DENSE || val > QPSolver::OSQP_SPARSE)
            return (MessageHandler::PrintMessage(INVALID_QPSOLVER,WARNING));

        qpSolver = (QPSolver) val;
        return SUCCESSFUL_RETURN;
    }


    void Options::setToDefault( ) {
        complementarityTolerance = 1.0e3 * Utilities::EPS;
        stationarityTolerance  = 1.0e6 * Utilities::EPS;
        initialPenaltyParameter = 0.01;
    	penaltyUpdateFactor  = 2.0;

        solveZeroPenaltyFirst = true;

        maxIterations = 1000;
        maxRho = 1e8;

        nDynamicPenalty = 3;
        etaDynamicPenalty = 0.9;

        printLevel = PrintLevel::INNER_LOOP_ITERATES;

        storeSteps = false;

        qpSolver = QPSolver::QPOASES_DENSE;
    }
}
