#include <pybind11/pybind11.h>

#include "OutputStatistics.hpp"


namespace LCQPow {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(OutputStatistics, m) {
  py::class_<OutputStatistics>(m, "OutputStatistics")
    .def(py::init<>())
    .def("getIterTotal", &OutputStatistics::getIterTotal)
    .def("getIterOuter", &OutputStatistics::getIterOuter)
    .def("getSubproblemIter", &OutputStatistics::getSubproblemIter)
    .def("getRhoOpt", &OutputStatistics::getRhoOpt)
    .def("getSolutionStatus", &OutputStatistics::getSolutionStatus)
    .def("getQPSolverExitFlag", &OutputStatistics::getQPSolverExitFlag)
    .def("getInnerIters", &OutputStatistics::getInnerIters)
    .def("getSubproblemIters", &OutputStatistics::getSubproblemIters)
    .def("getAccuSubproblemIters", &OutputStatistics::getAccuSubproblemIters)
    .def("getStepLength", &OutputStatistics::getStepLength)
    .def("getStepSize", &OutputStatistics::getStepSize)
    .def("getStatVals", &OutputStatistics::getStatVals)
    .def("getObjVals", &OutputStatistics::getObjVals)
    .def("getPhiVals", &OutputStatistics::getPhiVals)
    .def("getMeritVals", &OutputStatistics::getMeritVals);
}

} // namespace python
} // namespace LCQPow