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


#ifndef PLOTMANAGER_LCQPROBLEM_HPP
#define PLOTMANAGER_LCQPROBLEM_HPP

#include <matplotlibcpp.h>
#include <vector>

namespace lcqpOASES {

    enum LCQPNAME {
        WARM_UP = 0,
        IVOCP = 1,
        MOVING_MASSES = 2
    };

    class PlotManager {
        public:
            PlotManager(int _nV, int _nC, int _nComp, LCQPNAME lcqpname);

            ~PlotManager();

            void CreateIVOCPPlots(const double* const _xk, const double* const _lb, const double* const _ub);

            void CreateIVOCPTrajectoryPlot();

            void CreateIVOCPComplementarityPlot();

        private:
            // Variables describing the LCQP
            int nV, nC, nComp;

            // Debugging variables for IVOCP example
            double* xk;
            double* lb;
            double* ub;

            std::vector<int> indicesX;
            std::vector<int> indicesY;
            std::vector<int> indiceslam0;
    };
}


#endif  // PLOTMANAGER_LCQPROBLEM_HPP