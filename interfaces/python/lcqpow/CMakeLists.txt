
macro(pybind11_add_lcqpow_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_include_directories(
    ${MODULE} 
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${osqp_include}
  )
  target_link_libraries(
    ${MODULE} 
    PRIVATE
    PUBLIC ${PROJECT_NAME}-shared
    PRIVATE ${qpoases_lib} ${osqp_lib}
  )
endmacro()

pybind11_add_lcqpow_module(LCQProblem)
pybind11_add_lcqpow_module(Options)
pybind11_add_lcqpow_module(OutputStatistics)
pybind11_add_lcqpow_module(Utilities)