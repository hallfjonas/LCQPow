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

#include "PlotManager.hpp"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#include<vector>

namespace lcqpOASES {
    PlotManager::PlotManager(int _nV, int _nC, int _nComp, LCQPNAME lcqpname)
    {
        nV = _nV;
        nC = _nC;
        nComp = _nComp;

        xk = new double[nV]();

        switch (lcqpname) {
        case LCQPNAME::IVOCP:
            
            for (int i = 1; i < nV; i++) {
                if (i % 3 == 1)
                    indicesX.push_back(i);
                else if (i % 3 == 1)
                    indicesY.push_back(i);
                else
                    indiceslam0.push_back(i);
            }
            break;
        
        default:
            break;
        }
    }

    void PlotManager::CreateIVOCPPlots(const double* const _xk) 
    {
        memcpy(xk, _xk, (unsigned int)nV*sizeof(double));

        CreateIVOCPTrajectoryPlot();
        CreateIVOCPComplementarityPlot();
    }

    void PlotManager::CreateIVOCPTrajectoryPlot()
    {
        std::vector<double> xVals;
        for (unsigned int i = 0; i < indicesX.size(); i++) {
            xVals.push_back(xk[indicesX[i]]);
        }
        
        plt::plot(xVals);
        
        // Add graph title
        plt::title("Trajectory");

        plt::save("plots/trajectory.png");

    }

    void PlotManager::CreateIVOCPComplementarityPlot()
    {
        std::vector<double> xVals;
        std::vector<double> yVals;
        std::vector<double> lam0Vals;
        std::vector<double> lam1Vals;        

        for (unsigned int i = 0; i < indicesX.size(); i++) {
            xVals.push_back(xk[indicesX[i]]);
        }

        for (unsigned int i = 0; i < indicesY.size(); i++) {
            yVals.push_back(xk[indicesY[i]]);
        }

        for (unsigned int i = 0; i < indiceslam0.size(); i++) {
            lam0Vals.push_back(xk[indiceslam0[i]]);
            lam1Vals.push_back(xk[indicesX[i]] + xk[indiceslam0[i]]);
        }

        plt::named_plot("y", yVals);
        plt::named_plot("lam-", lam0Vals);
        plt::named_plot("lam+", lam1Vals);
        
        // Add graph title
        plt::title("Complementarity variables");

        // Enable legend.
        plt::legend();

        // Show image.
        plt::show();

        // Export
        plt::save("plots/complementarity_variables.png");
    }
}
