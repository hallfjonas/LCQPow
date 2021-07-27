% This file is part of LCQPow.
% LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
% Copyright (C) 2020 - 2021 by Jonas Hall et al.
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
%                min   1/2*x'Hx + x'g
%                s.t.  lbS1 <= S1*x <= ubS1
%                      lbS2 <= S2*x <= ubS2
%                         0  = x'*S1'*S2*x
%                        lb <=  x <= ub      {optional}
%                       lbA <= Ax <= ubA     {optional}
%
% I) Call
%
%    [x,{y, stats}] = LCQPow( H,g,S1,S2,lbS1,ubS1,lbS2,ubS2,{params} )
% or
%    [x,{y, stats}] = LCQPow( H,g,S1,S2,lbS1,ubS1,lbS2,ubS2,lb,ub,{params} )
% or
%    [x,{y, stats}] = LCQPow( H,g,S1,S2,lbS1,ubS1,lbS2,ubS2,A,lbA,ubA,{params} )
% or
%    [x,{y, stats}] = LCQPow( H,g,S1,S2,lbS1,ubS1,lbS2,ubS2,A,lbA,ubA,lb,ub,{params}).
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
%                  maxIterations : Maximum number of iterations to be performed.
%                     printLevel : The amount of output to be printed.
%
% III) The outputs consist of primal and dual solutions and a statistics struct:
%                              x : The primal solution (or last iterate on failed call)
%                              y : The dual solution (or last iterate on failed call)
%              stats.iters_total : Total number of iterations
%              stats.iters_outer : Total number of outer iterations
%         stats.iters_subproblem : Total number of iterations taken by the QP subsolver
%                  stats.rho_opt : Penalty parameter at solution (or its most recent value on failed call)
%             stats.elapsed_time : Elapsed solver time (without loading the data)
%                stats.exit_flag : Exit flag (0 on success, else some error according to the enum ReturnValue within Utilities.hpp)
%