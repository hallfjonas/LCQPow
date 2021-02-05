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

#include "qpOASES.hpp"
#include "Utilities.hpp"
#include <vector>

using qpOASES::QProblem;

/**
 *	\brief Implements the online active set strategy for QPs with general constraints.
 *
 *	A class for setting up and solving quadratic programs. The main feature is
 *	the possibily to use the newly developed online active set strategy for
 * 	parametric quadratic programming.
 *
 *	\author Jonas Hall
 *	\version 1.0
 *	\date 2020 - 2021
 */

namespace lcqpOASES {
	class LCQProblem
	{
		/*
		*	PUBLIC MEMBER FUNCTIONS
		*/
		public:
			/** Constructor which takes the QP dimension and Hessian type
			 *  information. If the Hessian is the zero (i.e. HST_ZERO) or the
			 *  identity matrix (i.e. HST_IDENTITY), respectively, no memory
			 *  is allocated for it and a NULL pointer can be passed for it
			 *  to the init() functions. */
			LCQProblem(	int _nV,	  							/**< Number of variables. */
						int _nC,		  						/**< Number of constraints. */
						int _nComp 								/**< Number of complementarity constraints. */
						);


			/** Initialises a QP problem with given QP data and tries to solve it
			 *	using at most nWSR iterations. Depending on the parameter constellation it: \n
			*	1. 0,    0,    0    : starts with xOpt = 0, yOpt = 0 and gB/gC empty (or all implicit equality bounds), \n
			*	2. xOpt, 0,    0    : starts with xOpt, yOpt = 0 and obtain gB/gC by "clipping", \n
			*	3. 0,    yOpt, 0    : starts with xOpt = 0, yOpt and obtain gB/gC from yOpt != 0, \n
			*	4. 0,    0,    gB/gC: starts with xOpt = 0, yOpt = 0 and gB/gC, \n
			*	5. xOpt, yOpt, 0    : starts with xOpt, yOpt and obtain gB/gC from yOpt != 0, \n
			*	6. xOpt, 0,    gB/gC: starts with xOpt, yOpt = 0 and gB/gC, \n
			*	7. xOpt, yOpt, gB/gC: starts with xOpt, yOpt and gB/gC (assume them to be consistent!)
			*
			*  Note: This function internally calls solveInitialLCQP for initialisation!
			*
			*	\return SUCCESSFUL_RETURN \n
						RET_INIT_FAILED \n
						RET_INIT_FAILED_CHOLESKY \n
						RET_INIT_FAILED_TQ \n
						RET_INIT_FAILED_HOTSTART \n
						RET_INIT_FAILED_INFEASIBILITY \n
						RET_INIT_FAILED_UNBOUNDEDNESS \n
						RET_MAX_NWSR_REACHED \n
						RET_INVALID_ARGUMENTS */
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

			/** Initialises a QP problem with given data to be read from files and solves it
			 *	using at most nWSR iterations. Depending on the parameter constellation it: \n
			*	1. 0,    0,    0    : starts with xOpt = 0, yOpt = 0 and gB/gC empty (or all implicit equality bounds), \n
			*	2. xOpt, 0,    0    : starts with xOpt, yOpt = 0 and obtain gB/gC by "clipping", \n
			*	3. 0,    yOpt, 0    : starts with xOpt = 0, yOpt and obtain gB/gC from yOpt != 0, \n
			*	4. 0,    0,    gB/gC: starts with xOpt = 0, yOpt = 0 and gB/gC, \n
			*	5. xOpt, yOpt, 0    : starts with xOpt, yOpt and obtain gB/gC from yOpt != 0, \n
			*	6. xOpt, 0,    gB/gC: starts with xOpt, yOpt = 0 and gB/gC, \n
			*	7. xOpt, yOpt, gB/gC: starts with xOpt, yOpt and gB/gC (assume them to be consistent!)
			*
			*  Note: This function internally calls solveInitialLCQP for initialisation!
			*
			*	\return SUCCESSFUL_RETURN \n
						RET_INIT_FAILED \n
						RET_INIT_FAILED_CHOLESKY \n
						RET_INIT_FAILED_TQ \n
						RET_INIT_FAILED_HOTSTART \n
						RET_INIT_FAILED_INFEASIBILITY \n
						RET_INIT_FAILED_UNBOUNDEDNESS \n
						RET_MAX_NWSR_REACHED \n
						RET_UNABLE_TO_READ_FILE \n
						RET_INVALID_ARGUMENTS */
			returnValue solve(	const char* const H_file,						/**< Name of file where Hessian matrix is stored. \n
																						If Hessian matrix is trivial, a NULL pointer can be passed. */
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

			/** Returns the dual solution vector of the LCQP (deep copy).
			 *	\return SUCCESSFUL_RETURN \n
						RET_QP_NOT_SOLVED */
			virtual algorithmStatus getPrimalSolution(	double* const xOpt	/**< Output: Primal solution vector (if QP has been solved). */
														) const;

			/** Returns the dual solution vector of the LCQP (deep copy).
			 *	\return SUCCESSFUL_RETURN \n
						RET_QP_NOT_SOLVED */
			virtual algorithmStatus getDualSolution(	double* const yOpt	/**< Output: Dual solution vector (if QP has been solved). */
														) const;

			/** Overrides current options with given ones.
			 *	\return SUCCESSFUL_RETURN */
			inline void setOptions(	const Options& _options	/**< New options. */
									);


		/*
		*	PROTECTED MEMBER FUNCTIONS
		*/
		protected:
			/** Frees all allocated memory.
			 *  \return SUCCESSFUL_RETURN */
			returnValue clear( );

			/** Copies all members from given rhs object.
			 *  \return SUCCESSFUL_RETURN */
			returnValue copy(	const LCQProblem& rhs	/**< Rhs object. */
								);


			/** Sets up dense internal QP data. If the current Hessian is trivial
			 *  (i.e. HST_ZERO or HST_IDENTITY) but a non-trivial one is given,
			 *  memory for Hessian is allocated and it is set to the given one.
			 *	\return SUCCESSFUL_RETURN \n
						RET_INVALID_ARGUMENTS \n
						RET_UNKNONW_BUG */
			returnValue setupLCQPdata(	const double* const _H, 		/**< Hessian matrix. \n
																			If Hessian matrix is trivial,a NULL pointer can be passed. */
										const double* const _g, 		/**< Gradient vector. */
										const double* const _A,  		/**< Constraint matrix. */
										const double* const _lb,		/**< Lower bound vector (on variables). \n
																			If no lower bounds exist, a NULL pointer can be passed. */
										const double* const _ub,		/**< Upper bound vector (on variables). \n
																			If no upper bounds exist, a NULL pointer can be passed. */
										const double* const _lbA,		/**< Lower constraints' bound vector. \n
																			If no lower constraints' bounds exist, a NULL pointer can be passed. */
										const double* const _ubA,		/**< Upper constraints' bound vector. \n
																			If no lower constraints' bounds exist, a NULL pointer can be passed. */
										const double* const _S1,      	/**< LHS of complementarity product. */
										const double* const _S2,      	/**< RHS of complementarity product. */
										const double* const _x0,		/**< Initial guess. */
										const double* const _y0			/**< Initial dual guess. */
										);

			/** Sets up internal QP data by loading it from files. If the current Hessian
			 *  is trivial (i.e. HST_ZERO or HST_IDENTITY) but a non-trivial one is given,
			 *  memory for Hessian is allocated and it is set to the given one.
			 *	\return SUCCESSFUL_RETURN \n
						RET_UNABLE_TO_OPEN_FILE \n
						RET_UNABLE_TO_READ_FILE \n
						RET_INVALID_ARGUMENTS \n
						RET_UNKNONW_BUG */
			returnValue setupLCQPdata(	const char* const H_file, 			/**< Name of file where Hessian matrix, of neighbouring QP to be solved, is stored. \n
																					If Hessian matrix is trivial,a NULL pointer can be passed. */
										const char* const g_file, 			/**< Name of file where gradient, of neighbouring QP to be solved, is stored. */
										const char* const A_file,			/**< Name of file where constraint matrix, of neighbouring QP to be solved, is stored. */
										const char* const lb_file, 			/**< Name of file where lower bounds, of neighbouring QP to be solved, is stored. \n
																					If no lower bounds exist, a NULL pointer can be passed. */
										const char* const ub_file, 			/**< Name of file where upper bounds, of neighbouring QP to be solved, is stored. \n
																					If no upper bounds exist, a NULL pointer can be passed. */
										const char* const lbA_file, 		/**< Name of file where lower constraints' bounds, of neighbouring QP to be solved, is stored. \n
																					If no lower constraints' bounds exist, a NULL pointer can be passed. */
										const char* const ubA_file,			/**< Name of file where upper constraints' bounds, of neighbouring QP to be solved, is stored. \n
																					If no upper constraints' bounds exist, a NULL pointer can be passed. */
										const char* const S1_file,			/**< Name of file where LHS of complementarity product is stored. */
										const char* const S2_file,			/**< Name of file where RHS of complementarity product is stored. */
										const char* const x0_file,			/**< Name of file where initial guess is stored. */
										const char* const y0_file			/**< Name of file where initial dual guess is stored. */
										);

			/** Prints concise information on the current iteration. */
			void printIteration( );

			/** Print header every once in a while. */
			void printHeader();

			/** Print line (mainly for printing header). */
			void printLine();


			/** Sets dense Hessian matrix of the QP.
			 *  If a null pointer is passed and
			 *  a) hessianType is HST_IDENTITY, nothing is done,
			 *  b) hessianType is not HST_IDENTITY, Hessian matrix is set to zero.
			 *	\return SUCCESSFUL_RETURN */
			inline returnValue setH(	const double* const H_new	/**< New dense Hessian matrix (with correct dimension!), a shallow copy is made. */
										);

										/** Changes gradient vector of the QP.
			 *	\return SUCCESSFUL_RETURN \n
			*			RET_INVALID_ARGUMENTS */
			inline returnValue setG(	const double* const g_new	/**< New gradient vector (with correct dimension!). */
										);

			/** Changes lower bound vector of the QP.
			 *	\return SUCCESSFUL_RETURN \n
			*			RET_QPOBJECT_NOT_SETUP */
			inline returnValue setLB(	const double* const lb_new	/**< New lower bound vector (with correct dimension!). */
										);

			/** Changes single entry of lower bound vector of the QP.
			 *	\return SUCCESSFUL_RETURN \n
			*			RET_QPOBJECT_NOT_SETUP \n
			*			RET_INDEX_OUT_OF_BOUNDS */
			inline returnValue setLB(	int number,	/**< Number of entry to be changed. */
										double value	/**< New value for entry of lower bound vector. */
										);

			/** Changes upper bound vector of the QP.
			 *	\return SUCCESSFUL_RETURN \n
			*			RET_QPOBJECT_NOT_SETUP */
			inline returnValue setUB(	const double* const ub_new	/**< New upper bound vector (with correct dimension!). */
										);

			/** Changes single entry of upper bound vector of the QP.
			 *	\return SUCCESSFUL_RETURN \n
			*			RET_QPOBJECT_NOT_SETUP \n
			*			RET_INDEX_OUT_OF_BOUNDS */
			inline returnValue setUB(	int number,	/**< Number of entry to be changed. */
										double value	/**< New value for entry of upper bound vector. */
										);

			inline returnValue setConstraints(  const double* const A_new,		/**< New constraint matrix. */
												const double* const S1_new,		/**< New lhs complementarity matrix. */
												const double* const S2_new,		/**< New rhs complementarity matrix. */
												const double* const lbA,		/**< New lower bounds for A. */
												const double* const ubA			/**< New upper bounds for A. */
												);

			/** Sets complementarity matrix (requires S1 and S2 to be set). \n
			 *	\return SUCCESSFUL_RETURN \n
			*			RET_INVALID_ARGUMENTS */
			returnValue setC( );

			/** Sets the initial guess x0 and y0. */
			inline returnValue setInitialGuess( const double* const _x0, const double* const _y0 );

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
			double* xnew;							/**< Current qpSubproblem solution. */
			double* pk;								/**< xnew - xk. */

			double alphak; 							/**< Optimal step length. */

			double* Qk;								/**< H + rho*C, required for stationarity and optimal step length. */
			double* statk;							/**< Stationarity of current iterate. */

			int outerIter;							/**< Outer iterate. */
			int innerIter;							/**< Inner iterate- */

			algorithmStatus algoStat;				/**< Status of algorithm. */

			QProblem qp;							/**< QP subproblem. */
			qpOASES::HessianType hessianType; 		/**< Hessian type (get this from initial solve). */
			int qpIterk;							/**< Iterations taken by qpSolver to solve subproblem. */
	};
}


#include "LCQProblem.ipp"

#endif	/* LCQPOASES_LCQPROBLEM_HPP */


/*
 *	end of file
 */
