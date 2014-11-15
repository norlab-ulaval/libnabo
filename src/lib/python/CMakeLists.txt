remake_python()
remake_project_get(PYTHON_MODULE_DESTINATION)
remake_add_library(
  pynabo PREFIX OFF
  *.cpp
  LINK nabo ${PYTHON_LIBRARY} ${Boost_LIBRARIES}
  COMPONENT python
  INSTALL ${PYTHON_MODULE_DESTINATION}/pynabo
)

remake_find_executable(python)
remake_test(Python pynabo
  PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/python/lib
    ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/bin/test/test.py
  DEPENDS pynabo
)
