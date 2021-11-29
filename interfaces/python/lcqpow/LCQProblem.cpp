#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "Eigen/Core"

#include "LCQProblem.hpp"


namespace LCQPow {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(LCQProblem, m) {
  py::class_<LCQProblem>(m, "LCQProblem")
    .def(py::init<>())
    .def(py::init<int, int, int>(), 
         py::arg("num_opt_vars"), py::arg("num_lin_constr"), 
         py::arg("num_compl_pairs"))
    .def("loadLCQP", [](LCQProblem& self, const Eigen::MatrixXd& H, 
                        const Eigen::VectorXd& g,  const Eigen::MatrixXd& S1, 
                        const Eigen::MatrixXd& S2, const Eigen::VectorXd& lbS1, 
                        const Eigen::VectorXd& ubS1, const Eigen::VectorXd& lbS2, 
                        const Eigen::VectorXd& ubS2, const Eigen::MatrixXd& A, 
                        const Eigen::VectorXd& lbA, const Eigen::VectorXd& ubA,
                        const Eigen::VectorXd& lb, const Eigen::VectorXd& ub,  
                        const Eigen::VectorXd& x0, const Eigen::VectorXd& y0) {
            return self.loadLCQP(H.data(), g.data(), S1.data(), S2.data(), 
                                lbS1.data(), ubS1.data(), lbS2.data(), ubS2.data(),
                                A.data(), lbA.data(), ubA.data(), 
                                lb.data(), ub.data(), x0.data(), y0.data());
          },
          py::arg("H"), py::arg("g"), py::arg("S1"), py::arg("S2"), 
          py::arg("lbS1"), py::arg("ubS1"), py::arg("lbS2"), py::arg("ubS2"), 
          py::arg("A"), py::arg("lbA"), py::arg("ubA"), 
          py::arg("lb"), py::arg("ub"), py::arg("x0"), py::arg("y0"))
    .def("loadLCQPFromFile", static_cast<ReturnValue (LCQProblem::*)(
          const char* const, const char* const, const char* const, 
          const char* const, const char* const, const char* const, 
          const char* const, const char* const, const char* const,
          const char* const, const char* const, const char* const,
          const char* const, const char* const, 
          const char* const)>(&LCQProblem::loadLCQP), 
          py::arg("H_file"), py::arg("g_file"), 
          py::arg("S1_file"), py::arg("S2_file"), 
          py::arg("lbS1_file")=nullptr, py::arg("ubS1_file")=nullptr, 
          py::arg("lbS2_file")=nullptr, py::arg("ubS2_file")=nullptr, 
          py::arg("A_file")=nullptr, py::arg("lbA_file")=nullptr, py::arg("ubA_file")=nullptr, 
          py::arg("lb_file")=nullptr, py::arg("ub_file")=nullptr, 
          py::arg("x0_file")=nullptr, py::arg("y0_file")=nullptr) 
    .def("runSolver", &LCQProblem::runSolver)
    .def("getPrimalSolution", &LCQProblem::getPrimalSolution)
    .def("getDualSolution", &LCQProblem::getDualSolution)
    .def("getNumerOfDuals", &LCQProblem::getNumerOfDuals)
    .def("getOutputStatistics", &LCQProblem::getOutputStatistics)
    .def("setOptions", &LCQProblem::setOptions);
}

} // namespace python
} // namespace LCQPow