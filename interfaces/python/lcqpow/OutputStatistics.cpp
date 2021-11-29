#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

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
    .def("getInnerIters", &OutputStatistics::getInnerItersStdVec)
    .def("getSubproblemIters", &OutputStatistics::getSubproblemItersStdVec)
    .def("getAccuSubproblemIters", &OutputStatistics::getAccuSubproblemItersStdVec)
    .def("getStepLength", &OutputStatistics::getStepLengthStdVec)
    .def("getStepSize", &OutputStatistics::getStepSizeStdVec)
    .def("getStatVals", &OutputStatistics::getStatValsStdVec)
    .def("getObjVals", &OutputStatistics::getObjValsStdVec)
    .def("getPhiVals", &OutputStatistics::getPhiValsStdVec)
    .def("getMeritVals", &OutputStatistics::getMeritValsStdVec);
}

} // namespace python
} // namespace LCQPow