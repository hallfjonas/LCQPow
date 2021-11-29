#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "Eigen/Core"

#include "LCQProblem.hpp"


namespace LCQPow {
namespace python {

namespace py = pybind11;


inline const double* getRawPtrFromEigenVectorXd(const Eigen::VectorXd& vec) {
  if (vec.size() > 0) return vec.data();
  else return nullptr;
}


inline const double* getRawPtrFromEigenMatrixXd(const Eigen::MatrixXd& mat) {
  if (mat.cols() > 0 && mat.rows() > 0)  return mat.data();
  else return nullptr;
}


PYBIND11_MODULE(LCQProblem, m) {
  py::class_<LCQProblem>(m, "LCQProblem")
    .def(py::init<>())
    .def(py::init<int, int, int>(), 
         py::arg("nV"), py::arg("nC"), py::arg("nComp"))
    .def("loadLCQP", [](LCQProblem& self, const Eigen::MatrixXd& H, 
                        const Eigen::VectorXd& g,  const Eigen::MatrixXd& S1, 
                        const Eigen::MatrixXd& S2, const Eigen::VectorXd& lbS1, 
                        const Eigen::VectorXd& ubS1, const Eigen::VectorXd& lbS2, 
                        const Eigen::VectorXd& ubS2, const Eigen::MatrixXd& A, 
                        const Eigen::VectorXd& lbA, const Eigen::VectorXd& ubA,
                        const Eigen::VectorXd& lb, const Eigen::VectorXd& ub,  
                        const Eigen::VectorXd& x0, const Eigen::VectorXd& y0) {
            return self.loadLCQP(H.data(), g.data(), S1.data(), S2.data(), 
                                 getRawPtrFromEigenVectorXd(lbS1),
                                 getRawPtrFromEigenVectorXd(ubS1),
                                 getRawPtrFromEigenVectorXd(lbS2),
                                 getRawPtrFromEigenVectorXd(ubS2),
                                 getRawPtrFromEigenMatrixXd(A),
                                 getRawPtrFromEigenVectorXd(lbA),
                                 getRawPtrFromEigenVectorXd(ubA),
                                 getRawPtrFromEigenVectorXd(lb),
                                 getRawPtrFromEigenVectorXd(ub),
                                 getRawPtrFromEigenVectorXd(x0),
                                 getRawPtrFromEigenVectorXd(y0));
          },
          py::arg("H"), py::arg("g"), py::arg("S1"), py::arg("S2"), 
          py::arg("lbS1")=Eigen::VectorXd::Zero(0), 
          py::arg("ubS1")=Eigen::VectorXd::Zero(0), 
          py::arg("lbS2")=Eigen::VectorXd::Zero(0), 
          py::arg("ubS2")=Eigen::VectorXd::Zero(0), 
          py::arg("A")=Eigen::MatrixXd::Zero(0, 0), 
          py::arg("lbA")=Eigen::VectorXd::Zero(0), 
          py::arg("ubA")=Eigen::VectorXd::Zero(0), 
          py::arg("lb")=Eigen::VectorXd::Zero(0), 
          py::arg("ub")=Eigen::VectorXd::Zero(0), 
          py::arg("x0")=Eigen::VectorXd::Zero(0), 
          py::arg("y0")=Eigen::VectorXd::Zero(0))
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
          py::arg("A_file")=nullptr, 
          py::arg("lbA_file")=nullptr, py::arg("ubA_file")=nullptr, 
          py::arg("lb_file")=nullptr, py::arg("ub_file")=nullptr, 
          py::arg("x0_file")=nullptr, py::arg("y0_file")=nullptr) 
    .def("runSolver", &LCQProblem::runSolver)
    .def("getPrimalSolution", [](const LCQProblem& self) {
            Eigen::VectorXd xOpt(Eigen::VectorXd::Zero(self.getNumerOfPrimals()));
            self.getPrimalSolution(xOpt.data());
            return xOpt;
         })
    .def("getDualSolution", [](const LCQProblem& self) {
            Eigen::VectorXd yOpt(Eigen::VectorXd::Zero(self.getNumerOfDuals()));
            self.getDualSolution(yOpt.data());
            return yOpt;
         })
    .def("getNumerOfDuals", &LCQProblem::getNumerOfDuals)
    .def("getOutputStatistics", &LCQProblem::getOutputStatistics)
    .def("setOptions", &LCQProblem::setOptions);
}

} // namespace python
} // namespace LCQPow