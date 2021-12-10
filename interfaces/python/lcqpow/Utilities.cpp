#include <pybind11/pybind11.h>

#include "Utilities.hpp"


namespace LCQPow {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(Utilities, m) {
  py::enum_<ReturnValue>(m, "ReturnValue", py::arithmetic())
    // Special values
    .value("NOT_YET_IMPLEMENTED",  ReturnValue::NOT_YET_IMPLEMENTED)
    .value("SUCCESSFUL_RETURN",  ReturnValue::SUCCESSFUL_RETURN)
    // Invalid arguments
    .value("INVALID_ARGUMENT",  ReturnValue::INVALID_ARGUMENT)
    .value("INVALID_PENALTY_UPDATE_VALUE",  ReturnValue::INVALID_PENALTY_UPDATE_VALUE)
    .value("INVALID_COMPLEMENTARITY_TOLERANCE",  ReturnValue::INVALID_COMPLEMENTARITY_TOLERANCE)
    .value("INVALID_INITIAL_PENALTY_VALUE",  ReturnValue::INVALID_INITIAL_PENALTY_VALUE)
    .value("INVALID_MAX_ITERATIONS_VALUE",  ReturnValue::INVALID_MAX_ITERATIONS_VALUE)
    .value("INVALID_STATIONARITY_TOLERANCE",  ReturnValue::INVALID_STATIONARITY_TOLERANCE)
    .value("INVALID_NUMBER_OF_OPTIM_VARS",  ReturnValue::INVALID_NUMBER_OF_OPTIM_VARS)
    .value("INVALID_NUMBER_OF_COMP_VARS ",  ReturnValue::INVALID_NUMBER_OF_COMP_VARS)
    .value("INVALID_NUMBER_OF_CONSTRAINT_VARS",  ReturnValue::INVALID_NUMBER_OF_CONSTRAINT_VARS)
    .value("INVALID_QPSOLVER",  ReturnValue::INVALID_QPSOLVER)
    .value("INVALID_OSQP_BOX_CONSTRAINTS",  ReturnValue::INVALID_OSQP_BOX_CONSTRAINTS)
    .value("INVALID_TOTAL_ITER_COUNT",  ReturnValue::INVALID_TOTAL_ITER_COUNT)
    .value("INVALID_TOTAL_OUTER_ITER",  ReturnValue::INVALID_TOTAL_OUTER_ITER)
    .value("IVALID_SUBPROBLEM_ITER",  ReturnValue::IVALID_SUBPROBLEM_ITER)
    .value("INVALID_RHO_OPT",  ReturnValue::INVALID_RHO_OPT)
    .value("INVALID_PRINT_LEVEL_VALUE",  ReturnValue::INVALID_PRINT_LEVEL_VALUE)
    .value("INVALID_OBJECTIVE_LINEAR_TERM",  ReturnValue::INVALID_OBJECTIVE_LINEAR_TERM)
    .value("INVALID_CONSTRAINT_MATRIX",  ReturnValue::INVALID_CONSTRAINT_MATRIX)
    .value("INVALID_COMPLEMENTARITY_MATRIX",  ReturnValue::INVALID_COMPLEMENTARITY_MATRIX)
    .value("INVALID_ETA_VALUE",  ReturnValue::INVALID_ETA_VALUE)
    .value("OSQP_INITIAL_PRIMAL_GUESS_FAILED",  ReturnValue::OSQP_INITIAL_PRIMAL_GUESS_FAILED)
    .value("OSQP_INITIAL_DUAL_GUESS_FAILED",  ReturnValue::OSQP_INITIAL_DUAL_GUESS_FAILED)
    .value("INVALID_LOWER_COMPLEMENTARITY_BOUND",  ReturnValue::INVALID_LOWER_COMPLEMENTARITY_BOUND)
    .value("INVALID_MAX_RHO_VALUE",  ReturnValue::INVALID_MAX_RHO_VALUE)
    // Algorithmic errors
    .value("MAX_ITERATIONS_REACHED",  ReturnValue::MAX_ITERATIONS_REACHED)
    .value("MAX_PENALTY_REACHED",  ReturnValue::MAX_PENALTY_REACHED)
    .value("INITIAL_SUBPROBLEM_FAILED",  ReturnValue::INITIAL_SUBPROBLEM_FAILED)
    .value("SUBPROBLEM_SOLVER_ERROR",  ReturnValue::SUBPROBLEM_SOLVER_ERROR)
    .value("FAILED_SYM_COMPLEMENTARITY_MATRIX",  ReturnValue::FAILED_SYM_COMPLEMENTARITY_MATRIX)
    .value("FAILED_SWITCH_TO_SPARSE",  ReturnValue::FAILED_SWITCH_TO_SPARSE)
    .value("FAILED_SWITCH_TO_DENSE",  ReturnValue::FAILED_SWITCH_TO_DENSE)
    .value("OSQP_WORKSPACE_NOT_SET_UP",  ReturnValue::OSQP_WORKSPACE_NOT_SET_UP)
    // Generic errors
    .value("LCQPOBJECT_NOT_SETUP",  ReturnValue::LCQPOBJECT_NOT_SETUP)
    .value("INDEX_OUT_OF_BOUNDS ",  ReturnValue::INDEX_OUT_OF_BOUNDS)
    .value("UNABLE_TO_READ_FILE ",  ReturnValue::UNABLE_TO_READ_FILE)
    // Sparse matrices
    .value("INVALID_INDEX_POINTER ",  ReturnValue::INVALID_INDEX_POINTER)
    .value("INVALID_INDEX_ARRAY ",  ReturnValue::INVALID_INDEX_ARRAY)
    .export_values();

  py::enum_<AlgorithmStatus>(m, "AlgorithmStatus", py::arithmetic())
    .value("PROBLEM_NOT_SOLVED", AlgorithmStatus::PROBLEM_NOT_SOLVED)
    .value("W_STATIONARY_SOLUTION", AlgorithmStatus::W_STATIONARY_SOLUTION)
    .value("C_STATIONARY_SOLUTION", AlgorithmStatus::C_STATIONARY_SOLUTION)
    .value("M_STATIONARY_SOLUTION", AlgorithmStatus::M_STATIONARY_SOLUTION)
    .value("S_STATIONARY_SOLUTION", AlgorithmStatus::S_STATIONARY_SOLUTION)
    .export_values();

  py::enum_<PrintLevel>(m, "PrintLevel", py::arithmetic())
    .value("NONE", PrintLevel::NONE)
    .value("OUTER_LOOP_ITERATES", PrintLevel::OUTER_LOOP_ITERATES)
    .value("INNER_LOOP_ITERATES", PrintLevel::INNER_LOOP_ITERATES)
    .value("SUBPROBLEM_SOLVER_ITERATES", PrintLevel::SUBPROBLEM_SOLVER_ITERATES)
    .export_values();

  py::enum_<QPSolver>(m, "QPSolver", py::arithmetic())
    .value("QPOASES_DENSE", QPSolver::QPOASES_DENSE)
    .value("QPOASES_SPARSE", QPSolver::QPOASES_SPARSE)
    .value("OSQP_SPARSE", QPSolver::OSQP_SPARSE)
    .export_values();
}

} // namespace python
} // namespace LCQPow