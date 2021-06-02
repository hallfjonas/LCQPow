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


#ifndef LCQPanther_MESSAGEHANDLER_HPP
#define LCQPanther_MESSAGEHANDLER_HPP

#include "Utilities.hpp"

namespace LCQPanther {

    class MessageHandler {

        public:
            /** Print a message.
             *
             * @param ret A return value which will impact the print message.
             *
             * @returns Simply passes the return value that was given. */
            static ReturnValue PrintMessage( ReturnValue ret );


            /** Print a solution status.
             *
             * @param algoStat The status to be printed.
             *
             * @returns Simply passes the algorithm status that was given. */
            static AlgorithmStatus PrintSolution( AlgorithmStatus algoStat );


            /** Prints a solution tag. */
            static void PrintSolutionLine( );
    };
}

#endif  // LCQPanther_MESSAGEHANDLER_HPP
