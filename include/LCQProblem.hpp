/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpOASES/LCQProblem.hpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Declaration of the LCQProblem class which is able to use the newly
 *	developed online active set strategy for parametric quadratic programming.
 */



#ifndef QPOASES_LCQPROBLEM_HPP
#define QPOASES_LCQPROBLEM_HPP

#include <qpOASES.hpp>


BEGIN_NAMESPACE_QPOASES


/**
 *	\brief Implements the online active set strategy for QPs with general constraints.
 *
 *	A class for setting up and solving quadratic programs. The main feature is
 *	the possibily to use the newly developed online active set strategy for
 * 	parametric quadratic programming.
 *
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 */
class LCQProblem
{
	/* allow SolutionAnalysis class to access private members */
	friend class SolutionAnalysis;

	/*
	 *	PUBLIC MEMBER FUNCTIONS
	 */
	public:
		/** Constructor which takes the QP dimension and Hessian type
		 *  information. If the Hessian is the zero (i.e. HST_ZERO) or the
		 *  identity matrix (i.e. HST_IDENTITY), respectively, no memory
		 *  is allocated for it and a NULL pointer can be passed for it
		 *  to the init() functions. */
		LCQProblem(	int_t _nV,	  							/**< Number of variables. */
					int_t _nC,		  						/**< Number of constraints. */
					int_t _nComp, 							/**< Number of complementarity constraints. */
					HessianType _hessianType = HST_UNKNOWN,	/**< Type of Hessian matrix. */
					BooleanType allocDenseMats = BT_TRUE	/**< Enable allocation of dense matrices. */
					);

		/** TODO: WRITE DESCRIPTION */
		returnValue solve(	SymmetricMatrix *_H,							/**< Hessian matrix (a shallow copy is made). */
							const real_t* const _g, 						/**< Gradient vector. */
							Matrix *_A,  									/**< Constraint matrix (a shallow copy is made). */
							const real_t* const _lb,						/**< Lower bound vector (on variables). \n
																					If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const _ub,						/**< Upper bound vector (on variables). \n
																					If no upper bounds exist, a NULL pointer can be passed. */
							const real_t* const _lbA,						/**< Lower constraints' bound vector. \n
																					If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const real_t* const _ubA,						/**< Upper constraints' bound vector. \n
																					If no lower constraints' bounds exist, a NULL pointer can be passed. */
							Matrix *_S1,          		         			/**< LHS of complementarity product. */
							Matrix *_S2,          		         			/**< RHS of complementarity product. */
							int_t& nWSR,									/**< Input: Maximum number of working set recalculations when using initial homotopy.
																					Output: Number of performed working set recalculations. */
							real_t* const cputime = 0,						/**< Input: Maximum CPU time allowed for QP initialisation. \n
																					Output: CPU time spent for QP initialisation (if pointer passed). */
							const real_t* const xOpt = 0,					/**< Optimal primal solution vector. \n
																					(If a null pointer is passed, the old primal solution is kept!) */
							const real_t* const yOpt = 0,					/**< Optimal dual solution vector. \n
																					(If a null pointer is passed, the old dual solution is kept!) */
							const Bounds* const guessedBounds = 0,			/**< Optimal working set of bounds for solution (xOpt,yOpt). \n
																					(If a null pointer is passed, all bounds are assumed inactive!) */
							const Constraints* const guessedConstraints = 0,/**< Optimal working set of constraints for solution (xOpt,yOpt). \n
																					(If a null pointer is passed, all constraints are assumed inactive!) */
							const real_t* const _R = 0						/**< Pre-computed (upper triangular) Cholesky factor of Hessian matrix.
																					The Cholesky factor must be stored in a real_t array of size nV*nV
																					in row-major format. Note: Only used if xOpt/yOpt and gB are NULL! \n
																					(If a null pointer is passed, Cholesky decomposition is computed internally!) */
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
		returnValue solve(	const real_t* const _H,							/**< Hessian matrix (a shallow copy is made). \n
																					 If Hessian matrix is trivial, a NULL pointer can be passed. */
							const real_t* const _g,							/**< Gradient vector. */
							const real_t* const _A,							/**< Constraint matrix (a shallow copy is made). */
							const real_t* const _lb,						/**< Lower bound vector (on variables). \n
																					 If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const _ub,						/**< Upper bound vector (on variables). \n
																					 If no upper bounds exist, a NULL pointer can be passed. */
							const real_t* const _lbA,						/**< Lower constraints' bound vector. \n
																					 If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const real_t* const _ubA,						/**< Upper constraints' bound vector. \n
																				 	If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const real_t* const _S1,               			/**< LHS of complementarity product. */
							const real_t* const _S2,               			/**< RHS of complementarity product. */
							int_t& nWSR,									/**< Input: Maximum number of working set recalculations when using initial homotopy.
																					Output: Number of performed working set recalculations. */
							real_t* const cputime = 0,						/**< Input: Maximum CPU time allowed for QP initialisation. \n
																					Output: CPU time spent for QP initialisation (if pointer passed). */
							const real_t* const xOpt = 0,					/**< Optimal primal solution vector. \n
																					(If a null pointer is passed, the old primal solution is kept!) */
							const real_t* const yOpt = 0,					/**< Optimal dual solution vector. \n
																					(If a null pointer is passed, the old dual solution is kept!) */
							const Bounds* const guessedBounds = 0,			/**< Optimal working set of bounds for solution (xOpt,yOpt). \n
																					(If a null pointer is passed, all bounds are assumed inactive!) */
							const Constraints* const guessedConstraints = 0,/**< Optimal working set of constraints for solution (xOpt,yOpt). \n
																					(If a null pointer is passed, all constraints are assumed inactive!) */
							const real_t* const _R = 0						/**< Pre-computed (upper triangular) Cholesky factor of Hessian matrix.
																					The Cholesky factor must be stored in a real_t array of size nV*nV
																					in row-major format. Note: Only used if xOpt/yOpt and gB are NULL! \n
																					(If a null pointer is passed, Cholesky decomposition is computed internally!) */
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
							const char* const A_file,						/**< Name of file where constraint matrix is stored. */
							const char* const lb_file,						/**< Name of file where lower bound vector. \n
																					If no lower bounds exist, a NULL pointer can be passed. */
							const char* const ub_file,						/**< Name of file where upper bound vector. \n
																					If no upper bounds exist, a NULL pointer can be passed. */
							const char* const lbA_file,						/**< Name of file where lower constraints' bound vector. \n
																					If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const char* const ubA_file,						/**< Name of file where upper constraints' bound vector. \n
																					If no upper constraints' bounds exist, a NULL pointer can be passed. */
							const char* const _S1,               			/**< Name of file where LHS of complementarity product is stored. */
							const char* const _S2,               			/**< Name of file where RHS of complementarity product is stored. */
							int_t& nWSR,									/**< Input: Maximum number of working set recalculations when using initial homotopy.
																					Output: Number of performed working set recalculations. */
							real_t* const cputime = 0,						/**< Input: Maximum CPU time allowed for QP initialisation. \n
																					Output: CPU time spent for QP initialisation (if pointer passed). */
							const real_t* const xOpt = 0,					/**< Optimal primal solution vector. \n
																					(If a null pointer is passed, the old primal solution is kept!) */
							const real_t* const yOpt = 0,					/**< Optimal dual solution vector. \n
																					(If a null pointer is passed, the old dual solution is kept!) */
							const Bounds* const guessedBounds = 0,			/**< Optimal working set of bounds for solution (xOpt,yOpt). \n
																					(If a null pointer is passed, all bounds are assumed inactive!) */
							const Constraints* const guessedConstraints = 0,/**< Optimal working set of constraints for solution (xOpt,yOpt). \n
																					(If a null pointer is passed, all constraints are assumed inactive!) */
							const char* const R_file = 0					/**< Name of the file where a pre-computed (upper triangular) Cholesky factor
																					of the Hessian matrix is stored. \n
																					(If a null pointer is passed, Cholesky decomposition is computed internally!) */
							);

		/** Returns the dual solution vector of the LCQP (deep copy).
		 *	\return SUCCESSFUL_RETURN \n
					RET_QP_NOT_SOLVED */
		virtual returnValue getPrimalSolution(	real_t* const xOpt	/**< Output: Primal solution vector (if QP has been solved). */
												) const;

		/** Returns the dual solution vector of the LCQP (deep copy).
		 *	\return SUCCESSFUL_RETURN \n
					RET_QP_NOT_SOLVED */
		virtual returnValue getDualSolution(	real_t* const yOpt	/**< Output: Dual solution vector (if QP has been solved). */
												) const;

		/** Overrides current options with given ones.
 		 *	\return SUCCESSFUL_RETURN */
		inline returnValue setOptions(	const Options& _options	/**< New options. */
										);

		/** Returns the print level.
		 *	\return Print level. */
		inline PrintLevel getPrintLevel( ) const;

		/** Changes the print level.
 		 *	\return SUCCESSFUL_RETURN */
		returnValue setPrintLevel(	PrintLevel _printlevel	/**< New print level. */
									);

		/** Returns the optimal objective function value.
		 *	\return finite value: Optimal objective function value (QP was solved) \n
		 			+infinity:	  QP was not yet solved */
		real_t getObjVal( ) const;

		/** Returns the objective function value at an arbitrary point x.
		 *	\return Objective function value at point x */
		real_t getObjVal(	const real_t* const _x	/**< Point at which the objective function shall be evaluated. */
							) const;


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

		/** Sets up internal QP data.
		 *	\return SUCCESSFUL_RETURN \n
					RET_INVALID_ARGUMENTS \n
					RET_UNKNONW_BUG */
		returnValue setupLCQPdata(	SymmetricMatrix *_H, 		/**< Hessian matrix. \n
																     If Hessian matrix is trivial,a NULL pointer can be passed. */
									const real_t* const _g, 	/**< Gradient vector. */
									Matrix *_A, 			 	/**< Constraint matrix. */
									const real_t* const _lb,	/**< Lower bound vector (on variables). \n
																	 If no lower bounds exist, a NULL pointer can be passed. */
									const real_t* const _ub,	/**< Upper bound vector (on variables). \n
																	 If no upper bounds exist, a NULL pointer can be passed. */
									const real_t* const _lbA,	/**< Lower constraints' bound vector. \n
																	 If no lower constraints' bounds exist, a NULL pointer can be passed. */
									const real_t* const _ubA,	/**< Upper constraints' bound vector. \n
																	 If no lower constraints' bounds exist, a NULL pointer can be passed. */
                                    Matrix *_S1,		      	/**< LHS of complementarity product. */
									Matrix *_S2 		   		/**< RHS of complementarity product. */
									);


		/** Sets up dense internal QP data. If the current Hessian is trivial
		 *  (i.e. HST_ZERO or HST_IDENTITY) but a non-trivial one is given,
		 *  memory for Hessian is allocated and it is set to the given one.
		 *	\return SUCCESSFUL_RETURN \n
					RET_INVALID_ARGUMENTS \n
					RET_UNKNONW_BUG */
		returnValue setupLCQPdata(	const real_t* const _H, 		/**< Hessian matrix. \n
																	     If Hessian matrix is trivial,a NULL pointer can be passed. */
									const real_t* const _g, 		/**< Gradient vector. */
									const real_t* const _A,  		/**< Constraint matrix. */
									const real_t* const _lb,		/**< Lower bound vector (on variables). \n
																		 If no lower bounds exist, a NULL pointer can be passed. */
									const real_t* const _ub,		/**< Upper bound vector (on variables). \n
																		 If no upper bounds exist, a NULL pointer can be passed. */
									const real_t* const _lbA,		/**< Lower constraints' bound vector. \n
																		 If no lower constraints' bounds exist, a NULL pointer can be passed. */
									const real_t* const _ubA,		/**< Upper constraints' bound vector. \n
																		 If no lower constraints' bounds exist, a NULL pointer can be passed. */
                                    const real_t* const _S1,      	/**< LHS of complementarity product. */
									const real_t* const _S2      	/**< RHS of complementarity product. */
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
									const char* const S2_file
									);

		/** Prints concise information on the current iteration.
		 *	\return  SUCCESSFUL_RETURN \n */
		returnValue printIteration(	int_t iter,							/**< Number of current iteration. */
									int_t BC_idx, 						/**< Index of blocking constraint. */
									SubjectToStatus BC_status,			/**< Status of blocking constraint. */
									BooleanType BC_isBound,				/**< Indicates if blocking constraint is a bound. */
									real_t homotopyLength,				/**< Current homotopy distance. */
									real_t penVal,						/**< Current penalty vaue. */
									BooleanType isFirstCall = BT_TRUE	/**< Indicating whether this is the first call for current QP. */
		 							);



		/** Sets Hessian matrix of the QP.
		 *	\return SUCCESSFUL_RETURN */
		inline returnValue setH(	SymmetricMatrix* H_new	/**< New Hessian matrix (a shallow copy is made). */
									);

		/** Sets dense Hessian matrix of the QP.
		 *  If a null pointer is passed and
		 *  a) hessianType is HST_IDENTITY, nothing is done,
		 *  b) hessianType is not HST_IDENTITY, Hessian matrix is set to zero.
		 *	\return SUCCESSFUL_RETURN */
		inline returnValue setH(	const real_t* const H_new	/**< New dense Hessian matrix (with correct dimension!), a shallow copy is made. */
									);

									/** Changes gradient vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		inline returnValue setG(	const real_t* const g_new	/**< New gradient vector (with correct dimension!). */
									);

		/** Changes lower bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP */
		inline returnValue setLB(	const real_t* const lb_new	/**< New lower bound vector (with correct dimension!). */
									);

		/** Changes single entry of lower bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP \n
		 *			RET_INDEX_OUT_OF_BOUNDS */
		inline returnValue setLB(	int_t number,	/**< Number of entry to be changed. */
									real_t value	/**< New value for entry of lower bound vector. */
									);

		/** Changes upper bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP */
		inline returnValue setUB(	const real_t* const ub_new	/**< New upper bound vector (with correct dimension!). */
									);

		/** Changes single entry of upper bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP \n
		 *			RET_INDEX_OUT_OF_BOUNDS */
		inline returnValue setUB(	int_t number,	/**< Number of entry to be changed. */
									real_t value	/**< New value for entry of upper bound vector. */
									);

		/** Sets constraint matrix of the QP. \n
			Note: Also internal vector Ax is recomputed!
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		inline returnValue setA(	Matrix *A_new	/**< New constraint matrix (a shallow copy is made). */
									);

		/** Sets dense constraint matrix of the QP. \n
			Note: Also internal vector Ax is recomputed!
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		inline returnValue setA(	const real_t* const A_new,	/**< New dense constraint matrix (with correct dimension!), a shallow copy is made. */
									const int_t nRows
									);


		/** Sets constraints' lower bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP */
		inline returnValue setLBA(	const real_t* const lbA_new	/**< New constraints' lower bound vector (with correct dimension!). */
									);

		/** Changes single entry of lower constraints' bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP \n
		 *			RET_INDEX_OUT_OF_BOUNDS */
		inline returnValue setLBA(	int_t number,	/**< Number of entry to be changed. */
									real_t value	/**< New value for entry of lower constraints' bound vector (with correct dimension!). */
									);

		/** Sets constraints' upper bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP */
		inline returnValue setUBA(	const real_t* const ubA_new	/**< New constraints' upper bound vector (with correct dimension!). */
									);

		/** Changes single entry of upper constraints' bound vector of the QP.
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_QPOBJECT_NOT_SETUP \n
		 *			RET_INDEX_OUT_OF_BOUNDS */
		inline returnValue setUBA(	int_t number,	/**< Number of entry to be changed. */
									real_t value	/**< New value for entry of upper constraints' bound vector (with correct dimension!). */
									);


		/** Returns number of complementarity constraints. \n
		 * 	\return SUCCESSFUL_RETURN  */
		inline int_t getNV( ) const;

				/** Returns number of complementarity constraints. \n
		 * 	\return SUCCESSFUL_RETURN  */
		inline int_t getNC( ) const;

		/** Returns number of complementarity constraints. \n
		 * 	\return SUCCESSFUL_RETURN  */
		inline int_t getNComp( ) const;


		/** Sets dense LHS complementarity matrix of the QP. \n
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		returnValue setComplementarities(	Matrix* S1_new,
											Matrix* S2_new
											);

		/** Sets dense LHS complementarity matrix of the QP. \n
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		returnValue setComplementarities(	const real_t* const S1_new,
											const real_t* const S2_new
											);

		/** Sets complementarity matrix (requires S1 and S2 to be set). \n
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		returnValue setC();


		/** Sets constraint matrix and vectors (requires S1 and S2 to be set if appendComplmenetarities is true). \n
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		returnValue setObjective(	SymmetricMatrix* H_new, 
									const real_t* const g_new
									);


		/** Sets constraint matrix and vectors (requires S1 and S2 to be set if appendComplmenetarities is true). \n
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		returnValue setConstraints(	const real_t* const A_new, 
									const real_t* const lbA_new,
									const real_t* const ubA_new
									);

		/** Sets constraint matrix and vectors (requires S1 and S2 to be set if appendComplmenetarities is true). \n
		 *	\return SUCCESSFUL_RETURN \n
		 *			RET_INVALID_ARGUMENTS */
		returnValue setConstraints(	Matrix* A_new, 
									const real_t* const lbA_new,
									const real_t* const ubA_new
									);

		/** Returns the NLP stationarity value. */
		real_t* getPenaltyStationarity(	const real_t* const x_eval,
										const real_t penVal,
										const real_t* const g_original
										);


	/*
	 *	PROTECTED MEMBER VARIABLES
	 */
	protected:
		BooleanType freeComplementarityMatrix;		/**< Flag indicating whether complementarity matrices should be de-allocated. */

		Matrix* S1;									/**< LHS of complementarity product. */
		Matrix* S2;									/**< RHS of complementarity product. */
		SymmetricMatrix* C;							/**< Complementarity matrix (S1'*S2 + S2'*S1). */

		Options options;							/**< Class for algorithmic options. */

	private:				
		int_t nV;									/**< Number of variables. */
		int_t nC;									/**< Number of constraints. */
		int_t nComp;								/**< Number of complementarity constraints. */

		BooleanType freeHessian;					/**< Flag indicating whether Hessian matrix should be de-allocated. */
		SymmetricMatrix* H;							/**< Objective Hessian term. */

		real_t* g;									/**< Objective linear term. */
		real_t* lb;									/**< Lower bound vector (on variables). */
		real_t* ub;									/**< Upper bound vector (on variables). */

		BooleanType freeConstraintMatrix;			/**< Flag indicating whether constraint matrix should be de-allocated. */
		Matrix* A;									/**< Constraint matrix. */
		real_t* lbA;								/**< Lower bound vector (on constraints). */
		real_t* ubA;								/**< Upper bound vector (on constraints). */

		QProblem qp;								/**< QP subproblem. */
};


END_NAMESPACE_QPOASES

#include <qpOASES/LCQProblem.ipp>

#endif	/* QPOASES_LCQPROBLEM_HPP */


/*
 *	end of file
 */
