/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
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


#include "MessageHandler.hpp"


namespace LCQPow {

    ReturnValue MessageHandler::PrintMessage( ReturnValue ret) {

        switch (ret) {
            case SUCCESSFUL_RETURN:
                break;

            case NOT_YET_IMPLEMENTED:
                printf("This method has not yet been implemented.\n");
                break;

            case LCQPOBJECT_NOT_SETUP:
                printf("ERROR: The LCQP object has not been set up correctly.\n");
                break;

            case INDEX_OUT_OF_BOUNDS:
                printf("ERROR: Index out of bounds.\n");
                break;

            case SUBPROBLEM_SOLVER_ERROR:
                printf("ERROR: The subproblem solver produced an error.\n");
                break;

            case UNABLE_TO_READ_FILE:
                printf("ERROR: Unable to read file.\n");
                break;

            case MAX_ITERATIONS_REACHED:
                printf("ERROR: Maximum number of iterations reached.\n");
                break;

            case MAX_PENALTY_REACHED:
                printf("ERROR: Maxium penalty value reached.\n");
                break;

            case INITIAL_SUBPROBLEM_FAILED:
                printf("ERROR: Failed to solve initial QP.\n");
                break;

            case INVALID_ARGUMENT:
                printf("ERROR: Invalid argument passed.\n");
                break;

            case INVALID_NUMBER_OF_OPTIM_VARS:
                printf("ERROR: Invalid optimization variable dimension passed (required to be > 0).\n");
                break;

            case INVALID_NUMBER_OF_COMP_VARS:
                printf("ERROR: Invalid complementarity dimension passed (required to be > 0).\n");
                break;

            case INVALID_NUMBER_OF_CONSTRAINT_VARS:
                printf("ERROR: Invalid number of optimization variables passed (required to be >= 0).\n");
                break;

            case INVALID_QPSOLVER:
                printf("ERROR: Invalid QPSolver passed.\n");
                break;

            case INVALID_COMPLEMENTARITY_TOLERANCE:
                printf("WARNING: Ignoring invalid complementarity tolerance.\n");
                break;

            case INVALID_INITIAL_PENALTY_VALUE:
                printf("WARNING: Invalid argument passed (initial penalty value).\n");
                break;

            case INVALID_PENALTY_UPDATE_VALUE:
                printf("WARNING: Ignoring invalid penalty update value.\n");
                break;

            case INVALID_MAX_ITERATIONS_VALUE:
                printf("WARNING: Ignoring invalid number of maximum iterations.\n");
                break;

            case INVALID_MAX_RHO_VALUE:
                printf("WARNING: Ignoring invalid number of maximum penalty value.\n");
                break;

            case INVALID_STATIONARITY_TOLERANCE:
                printf("WARNING: Ignoring invalid stationarity tolerance.\n");
                break;

            case INVALID_INDEX_POINTER:
                printf("ERROR: Invalid index pointer passed in csc format.\n");
                break;

            case INVALID_INDEX_ARRAY:
                printf("ERROR: Invalid index array passed in csc format.\n");
                break;

            case INVALID_OSQP_BOX_CONSTRAINTS:
                printf("ERROR: Invalid constraints passed to OSQP solver: This solver does not handle box constraints, please pass them through linear constraints.\n");
                break;

            case INVALID_TOTAL_ITER_COUNT:
                printf("ERROR: Invalid total number of iterations delta passed to output statistics (must be non-negative integer).\n");
                break;

            case INVALID_TOTAL_OUTER_ITER:
                printf("ERROR: Invalid total number of outer iterations delta passed to output statistics (must be non-negative integer).\n");
                break;

            case IVALID_SUBPROBLEM_ITER:
                printf("ERROR: Invalid total number of subproblem solver iterates delta passed to output statistics (must be non-negative integer).\n");
                break;

            case INVALID_RHO_OPT:
                printf("ERROR: Invalid rho value at solution passed to output statistics (must be positive double).\n");
                break;

            case INVALID_PRINT_LEVEL_VALUE:
                printf("WARNING: Ignoring invalid integer to be parsed to print level passed (must be in range of enum).\n");
                break;

            case INVALID_OBJECTIVE_LINEAR_TERM:
                printf("ERROR: Invalid objective linear term passed (must be a double array of length n).\n");
                break;

            case INVALID_CONSTRAINT_MATRIX:
                printf("ERROR: Invalid constraint matrix passed (matrix was null pointer but number of constraints is positive).\n");
                break;

            case INVALID_COMPLEMENTARITY_MATRIX:
                printf("ERROR: Invalid complementarity matrix passed (can not be null pointer).\n");
                break;

            case INVALID_ETA_VALUE:
                printf("ERROR: Invalid etaComplHist value, which describes the fraction of loss required for complementarity progress (must be in (0,1)).");
                break;

            case OSQP_INITIAL_PRIMAL_GUESS_FAILED:
                printf("ERROR: OSQP failed to use the primal initial guess.\n");
                break;

            case OSQP_INITIAL_DUAL_GUESS_FAILED:
                printf("ERROR: OSQP failed to use the dual initial guess.\n");
                break;

            case INVALID_LOWER_COMPLEMENTARITY_BOUND:
                printf("ERROR: Lower complementarity bound must be bounded below.\n");
                break;

            case FAILED_SYM_COMPLEMENTARITY_MATRIX:
                printf("Failed to compute the symmetric complementarity matrix C.\n");
                break;

            case FAILED_SWITCH_TO_SPARSE:
                printf("Failed to switch to sparse mode (a to be created sparse matrix was nullpointer).\n");
                break;

            case FAILED_SWITCH_TO_DENSE:
                printf("Failed to switch to dense mode (an array to be created was nullpointer).\n");
                break;

            case OSQP_WORKSPACE_NOT_SET_UP:
                printf("OSQP Workspace is not set up (please check for OSQP errors).\n");
                break;
        }

        fflush(stdout);

        return ret;
    }


    AlgorithmStatus MessageHandler::PrintSolution( AlgorithmStatus algoStat ) {

        if ( algoStat != PROBLEM_NOT_SOLVED)
            printf("\n\n#################################\n");

        switch (algoStat) {
            case PROBLEM_NOT_SOLVED:
                printf("The LCQP has not been solved.\n");
                break;

            case W_STATIONARY_SOLUTION:
                printf("## W-Stationary solution found ##\n");
                break;

            case C_STATIONARY_SOLUTION:
                printf("## C-Stationary solution found ##\n");
                break;

            case M_STATIONARY_SOLUTION:
                printf("## M-Stationary solution found ##\n");
                break;

            case S_STATIONARY_SOLUTION:
                printf("## S-Stationary solution found ##\n");
                break;
        }

        if ( algoStat != PROBLEM_NOT_SOLVED)
            printf("#################################\n\n");

        return algoStat;
    }
}

