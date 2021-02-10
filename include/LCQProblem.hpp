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


#ifndef LCQPOASES_LCQPROBLEM_HPP
#define LCQPOASES_LCQPROBLEM_HPP

#include <Utilities.hpp>
#include <Subsolver.hpp>

#include <qpOASES.hpp>
#include <vector>

using qpOASES::QProblem;

namespace lcqpOASES {
	class LCQProblem
	{
		/*
		*	PUBLIC MEMBER FUNCTIONS
		*/
		public:
			/** Default constructor. */
			LCQProblem( );

			/** TODO: Write description. */
			LCQProblem(	int _nV,	  							/**< Number of variables. */
						int _nC,		  						/**< Number of constraints. */
						int _nComp 								/**< Number of complementarity constraints. */
						);


			/** TODO: Write description. */
			returnValue solve(	const double* const _H,							/**< Hessian matrix. */
								const double* const _g,							/**< Gradient vector. */
								const double* const _lb,						/**< Lower bound vector (on variables). \n
																						If no lower bounds exist, a NULL pointer can be passed. */
								const double* const _ub,						/**< Upper bound vector (on variables). \n
																						If no upper bounds exist, a NULL pointer can be passed. */
								const double* const _S1,               			/**< LHS of complementarity product. */
								const double* const _S2,               			/**< RHS of complementarity product. */
								const double* const _A = 0,						/**< Constraint matrix.
																						If no constraints exist, a NULL pointer can be passed. */
								const double* const _lbA = 0,					/**< Lower constraints' bound vector. \n
																						If no lower constraints' bounds exist, a NULL pointer can be passed. */
								const double* const _ubA = 0,					/**< Upper constraints' bound vector. \n
																						If no lower constraints' bounds exist, a NULL pointer can be passed. */
								const double* const xOpt = 0,					/**< Optimal primal solution vector. \n
																						(If a null pointer is passed, the old primal solution is kept!) */
								const double* const yOpt = 0					/**< Optimal dual solution vector. \n
																						(If a null pointer is passed, the old dual solution is kept!) */
								);

			/** TODO: Write description. */
			returnValue solve(	const char* const H_file,						/**< Name of file where Hessian matrix is stored. */
								const char* const g_file,						/**< Name of file where gradient vector is stored. */
								const char* const lb_file,						/**< Name of file where lower bound vector. \n
																						If no lower bounds exist, a NULL pointer can be passed. */
								const char* const ub_file,						/**< Name of file where upper bound vector. \n
																						If no upper bounds exist, a NULL pointer can be passed. */
								const char* const S1_file,             			/**< Name of file where LHS of complementarity product is stored. */
								const char* const S2_file,             			/**< Name of file where RHS of complementarity product is stored. */
								const char* const A_file = 0,					/**< Name of file where constraint matrix is stored.
																						If no constraints exist, a NULL pointer can be passed. */
								const char* const lbA_file = 0,					/**< Name of file where lower constraints' bound vector. \n
																						If no lower constraints' bounds exist, a NULL pointer can be passed. */
								const char* const ubA_file = 0,					/**< Name of file where upper constraints' bound vector. \n
																						If no upper constraints' bounds exist, a NULL pointer can be passed. */
								const char* const x0_file = 0,					/**< Optimal primal solution vector. \n
																						(If a null pointer is passed, the old primal solution is kept!) */
								const char* const y0_file = 0					/**< Optimal dual solution vector. \n
																						(If a null pointer is passed, the old dual solution is kept!) */
								);

			/** TODO: Write description. */
			returnValue solve(	double* _H_data,						/**< Non-zero Hessian matrix values. */
								int _H_nnx,								/**< Number of non-zero values of Hessian. */
								int* _H_i,								/**< Row indicies of non-zero values. */
								int* _H_p,								/**< Pointer to column starts. */
								double* _g,								/**< Gradient vector. */
								double* _lb,							/**< Lower bound vector (on variables). */
								double* _ub,							/**< Upper bound vector (on variables). */
								double* _S1_data,          				/**< Non-zero values of LHS of complementarity product. */
								int _S1_nnx,							/**< Number of non-zero values of S1. */
								int* _S1_i,								/**< Row indicies of non-zero values. */
								int* _S1_p,								/**< Pointer to column starts. */
								double* _S2_data,	          			/**< Non-zero values of RHS of complementarity product. */
								int _S2_nnx,							/**< Number of non-zero values of S2. */
								int* _S2_i,								/**< Row indicies of non-zero values. */
								int* _S2_p,								/**< Pointer to column starts. */
								double* _A_data = 0,					/**< Constraint matrix. */
								int _A_nnx = 0,							/**< Number of non-zero values of Hessian. */
								int* _A_i = 0,							/**< Row indicies of non-zero values. */
								int* _A_p = 0,							/**< Pointer to column starts. */
								double* _lbA = 0,						/**< Lower constraints' bound vector. */
								double* _ubA = 0,						/**< Upper constraints' bound vector. */
								double* _x0 = 0,						/**< Optimal primal solution vector. */
								double* _y0 = 0							/**< Optimal dual solution vector. */
								);

			/** TODO: Write description. */
			returnValue solve(	const char* const H_data_file,					/**< Non-zero Hessian matrix values. */
								const char* const H_i_file,						/**< Row indicies of non-zero values. */
								const char* const H_p_file,						/**< Pointer to column starts. */
								const char* const g_file,						/**< Gradient vector. */
								const char* const lb_file,						/**< Lower bound vector (on variables). */
								const char* const ub_file,						/**< Upper bound vector (on variables). */
								const char* const S1_data_file,          		/**< Non-zero values of LHS of complementarity product. */
								const char* const S1_i_file,					/**< Row indicies of non-zero values. */
								const char* const S1_p_file,					/**< Pointer to column starts. */
								const char* const S2_data_file,          		/**< Non-zero values of RHS of complementarity product. */
								const char* const S2_i_file,					/**< Row indicies of non-zero values. */
								const char* const S2_p_file,					/**< Pointer to column starts. */
								const char* const A_data_file = 0,				/**< Constraint matrix. */
								const char* const A_i_file = 0,					/**< Row indicies of non-zero values. */
								const char* const A_p_file = 0,					/**< Pointer to column starts. */
								const char* const lbA_file = 0,					/**< Lower constraints' bound vector. */
								const char* const ubA_file = 0,					/**< Upper constraints' bound vector. */
								const char* const x0_file = 0,					/**< Optimal primal solution vector. */
								const char* const y0_file = 0					/**< Optimal dual solution vector. */
								);

			/** TODO: Write description. */
			virtual algorithmStatus getPrimalSolution(	double* const xOpt	/**< Output: Primal solution vector (if QP has been solved). */
														) const;

			/** TODO: Write description. */
			virtual algorithmStatus getDualSolution(	double* const yOpt	/**< Output: Dual solution vector (if QP has been solved). */
														) const;

			/** TODO: Write description. */
			inline void setOptions(	const Options& _options	/**< New options. */
									);


		/*
		*	PROTECTED MEMBER FUNCTIONS
		*/
		protected:
			/** TODO: Write description. */
			returnValue clear( );


			/** TODO: Write description. */
			returnValue copy(	const LCQProblem& rhs	/**< Rhs object. */
								);

			/** Prints concise information on the current iteration. */
			void printIteration( );

			/** Print header every once in a while. */
			void printHeader();

			/** Print line (mainly for printing header). */
			void printLine();


			/** TODO: Write description. */
			inline returnValue setH(	const double* const H_new	/**< New dense Hessian matrix (with correct dimension!), a shallow copy is made. */
										);

			inline returnValue setH(	double* H_data, 
										int H_nnx,
										int* H_i,
										int* H_p 
										);

			/** TODO: Write description. */
			inline returnValue setG(	const double* const g_new	/**< New gradient vector (with correct dimension!). */
										);

			/** TODO: Write description. */
			inline returnValue setLB(	const double* const lb_new	/**< New lower bound vector (with correct dimension!). */
										);

			/** TODO: Write description. */
			inline returnValue setLB(	int number,	/**< Number of entry to be changed. */
										double value	/**< New value for entry of lower bound vector. */
										);

			/** TODO: Write description. */
			inline returnValue setUB(	const double* const ub_new	/**< New upper bound vector (with correct dimension!). */
										);

			/** TODO: Write description. */
			inline returnValue setUB(	int number,	/**< Number of entry to be changed. */
										double value	/**< New value for entry of upper bound vector. */
										);

			/** TODO: Write description. */
			returnValue setConstraints( const double* const S1_new,		/**< New lhs complementarity matrix. */
										const double* const S2_new,		/**< New rhs complementarity matrix. */
										const double* const A_new,		/**< New constraint matrix. */
										const double* const lbA,		/**< New lower bounds for A. */
										const double* const ubA			/**< New upper bounds for A. */
										);

			/** TODO: Write description. */
			returnValue setConstraints(	double* S1_data,		/**< New lhs complementarity matrix. */
										int S1_nnx,
										int* S1_i,
										int* S1_p,
										double* S2_data,		/**< New rhs complementarity matrix. */
										int S2_nnx,
										int* S2_i,
										int* S2_p,
										double* A_data,		/**< New constraint matrix. */
										int A_nnx,
										int* A_i,
										int* A_p,
										double* lbA,		/**< New lower bounds for A. */
										double* ubA			/**< New upper bounds for A. */
										);

			/** TODO: Write description. */
			returnValue setC( );

			/** Sets the initial guess x0 and y0. */
			inline returnValue setInitialGuess( const double* const _x0, const double* const _y0 );


			/** TODO: Write description. */
			inline returnValue setSparseMatrix(	const double* const _M_data,
												const int _M_nnx,
												const int* const _M_i,
												const int* const _M_p,
												qpOASES::SparseMatrix Mat
												);

			/** TODO: Write description. */
			inline returnValue setSparseMatrix(	const double* const _M_data,
												const int _M_nnx,
												const int* const _M_i,
												const int* const _M_p,
												qpOASES::SymSparseMat Mat
												);

			/** Returns the NLP stationarity value. */
			double* getPenaltyStationarity(	const double* const x_eval,
											const double penVal,
											const double* const g_original
											);


		/*
		*	PROTECTED MEMBER VARIABLES
		*/
		protected:
			Options options;						/**< Class for algorithmic options. */

		private:

			/** After problem is set up, call this function and solve LCQP. */
			returnValue runSolver( );

			/** Called in runSolver to initialize variables. */
			void initializeSolver( );

			/** Solves the qp subproblem wrt H, gk, A, S1, S2 . */
			returnValue solveQPSubproblem( bool initialSolve );

			/** Check outer stationarity at current iterate xk. */
			bool stationarityCheck( );

			/** Check satisfaction of complementarity value. */
			bool complementarityCheck( );

			/** Perform penalty update. */
			void updatePenalty( );

			/** Get optimal step length. */
			void getOptimalStepLength( );

			/** Update xk and gk. */
			void updateStep( );

			/** Gradient perturbation method. */
			void perturbGradient( );

			/** Transform the dual variables from penalty form to LCQP form. */
			void transformDuals( );

			/** Determine stationarity type of optimal solution. */
			void determineStationarityType( );

			/** Get indices of weak complementarites. */
			std::vector<int> getWeakComplementarities( );

			int nV;									/**< Number of variables. */
			int nC;									/**< Number of constraints. */
			int nComp;								/**< Number of complementarity constraints. */

			double* H;								/**< Objective Hessian term. */

			double* g;								/**< Objective linear term. */
			double* lb;								/**< Lower bound vector (on variables). */
			double* ub;								/**< Upper bound vector (on variables). */

			double* A;								/**< Constraint matrix. */
			double* lbA;							/**< Lower bound vector (on constraints). */
			double* ubA;							/**< Upper bound vector (on constraints). */

			double* S1;								/**< LHS of complementarity product. */
			double* S2;								/**< RHS of complementarity product. */
			double* C;								/**< Complementarity matrix (S1'*S2 + S2'*S1). */

			double rho; 							/**< Current penalty value. */

			double* gk;								/**< Current objective linear term. */

			double* xk;								/**< Current primal iterate. */
			double* yk;								/**< Current dual vector. */
			double* yk_A;							/**< Current dual vector w.r.t A. */
			double* xnew;							/**< Current qpSubproblem solution. */
			double* pk;								/**< xnew - xk. */

			double alphak; 							/**< Optimal step length. */

			double* Qk;								/**< H + rho*C, required for stationarity and optimal step length. */
			double* statk;							/**< Stationarity of current iterate. */

			int outerIter;							/**< Outer iterate. */
			int innerIter;							/**< Inner iterate- */

			int qpIterk;						/**< Iterations taken by qpSolver to solve subproblem. */

			algorithmStatus algoStat;				/**< Status of algorithm. */

			// Sparse matrices
			csc* H_sparse;
			csc* A_sparse;
			
			// Subproblem solvers
            Subsolver subsolver;

	};
}


#include "LCQProblem.ipp"

#endif	/* LCQPOASES_LCQPROBLEM_HPP */


/*
 *	end of file
 */
