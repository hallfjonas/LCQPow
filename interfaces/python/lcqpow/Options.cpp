#include <pybind11/pybind11.h>

#include "Options.hpp"


namespace LCQPow {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(Options, m) {
  py::class_<Options>(m, "Options")
    .def(py::init<>())
    .def(py::init<const Options&>(), py::arg("rhs"))
    .def("setToDefault", &Options::setToDefault)
    .def("getStationarityTolerance", &Options::getStationarityTolerance)
    .def("getComplementarityTolerance", &Options::getComplementarityTolerance)
    .def("setComplementarityTolerance", &Options::setComplementarityTolerance)
    .def("getInitialPenaltyParameter", &Options::getInitialPenaltyParameter)
    .def("setInitialPenaltyParameter", &Options::setInitialPenaltyParameter)
    .def("getPenaltyUpdateFactor", &Options::getPenaltyUpdateFactor)
    .def("setPenaltyUpdateFactor", &Options::setPenaltyUpdateFactor)
    .def("getSolveZeroPenaltyFirst", &Options::getSolveZeroPenaltyFirst)
    .def("setSolveZeroPenaltyFirst", &Options::setSolveZeroPenaltyFirst)
    .def("getMaxIterations", &Options::getMaxIterations)
    .def("setMaxIterations", &Options::setMaxIterations)
    .def("getMaxRho", &Options::getMaxRho)
    .def("setMaxRho", &Options::setMaxRho)
    .def("getNDynamicPenalty", &Options::getNDynamicPenalty)
    .def("setNDynamicPenalty", &Options::setNDynamicPenalty)
    .def("getEtaDynamicPenalty", &Options::getEtaDynamicPenalty)
    .def("setEtaDynamicPenalty", &Options::setEtaDynamicPenalty)
    .def("getPrintLevel", &Options::getPrintLevel)
    .def("setPrintLevel", static_cast<ReturnValue (Options::*)(PrintLevel)>(&Options::setPrintLevel))
    .def("setPrintLevel", static_cast<ReturnValue (Options::*)(int)>(&Options::setPrintLevel))
    .def("getStoreSteps", &Options::getStoreSteps)
    .def("getQPSolver", &Options::getQPSolver)
    .def("setQPSolver", static_cast<ReturnValue (Options::*)(QPSolver)>(&Options::setQPSolver))
    .def("setQPSolver", static_cast<ReturnValue (Options::*)(int)>(&Options::setQPSolver));
}

} // namespace python
} // namespace LCQPow