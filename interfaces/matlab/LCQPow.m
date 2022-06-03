% This file is part of LCQPow.
% LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
% Copyright (C) 2020 - 2022 by Jonas Hall et al.
%
% LCQPow is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
% LCQPow is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with LCQPow; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%---------------------------------------------------------------------------------
%
% LCQPow is intended for solving quadratic programs with
% linear complementarity constraints of the form
%
%           minimize   1/2*x'Qx + x'g
%             s.t.    0  = x'*L'*R*x
%                  lbL <= L*x <= ubL
%                  lbR <= R*x <= ubR
%                   lbA <=  Ax  <= ubA     {optional}
%                    lb <=   x  <= ub      {optional}
%
% I) Call
%
%    [x,{y, stats}] = LCQPow( Q,g,L,R,lbL,ubL,lbR,ubR,{params} )
% or
%    [x,{y, stats}] = LCQPow( Q,g,L,R,lbL,ubL,lbR,ubR,lb,ub,{params} )
% or
%    [x,{y, stats}] = LCQPow( Q,g,L,R,lbL,ubL,lbR,ubR,A,lbA,ubA,{params} )
% or
%    [x,{y, stats}] = LCQPow( Q,g,L,R,lbL,ubL,lbR,ubR,A,lbA,ubA,lb,ub,{params}).
%
% II) The optional params struct may contain the following fields:
%
%                             x0 : Initial primal guess.
%                             y0 : Initial dual guess.
%          stationarityTolerance : Tolerance for 1-Norm of stationarity violation.
%       complementarityTolerance : Complementarity tolerance.
%        initialPenaltyParameter : Start value for complementarity penalty term.
%            penaltyUpdateFactor : Factor for updating penaltised complementarity term.
%          solveZeroPenaltyFirst : Flag indicating whether first QP should ignore penalization.
%                    perturbStep : Flag indicating whether to perform step perturbation.
%                  maxIterations : Maximum number of iterations to be performed.
%                         maxRho : Maximum penalty value.
%                     printLevel : The amount of output to be printed.
%                     nDynamicPenalty : The number of complementarity values to be compared in Leyffer check.
%                   etaDynamicPenalty : Complementarity reduction factor required in at least one of the lastet nDynamicPenalty steps.
%
% III) The outputs consist of primal and dual solutions and a statistics struct:
%                              x : The primal solution (or last iterate on failed call)
%                              y : The dual solution (or last iterate on failed call)
%              stats.iters_total : Total number of iterations
%              stats.iters_outer : Total number of outer iterations
%         stats.iters_subproblem : Total number of iterations taken by the QP subsolver
%                  stats.rho_opt : Penalty parameter at solution (or its most recent value on failed call)
%             stats.elapsed_time : Elapsed solver time (without interface overhead)
%                stats.exit_flag : Exit flag (0 on success, else some error according to the enum ReturnValue within Utilities.hpp)
%             stats.qp_exit_flag : A flag indicating the most recent status flag of the QP solver.
%