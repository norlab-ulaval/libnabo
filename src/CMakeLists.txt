remake_find_package(
  Boost REQUIRED
  COMPONENTS python
  RESULT_VAR Boost_FOUND
)
remake_find_package(eigen3 CONFIG)
if(EIGEN3_FOUND)
  remake_project_set(EIGEN_VERSION 3 CACHE INTERNAL
    "Version of the eigen linear algebra library")
else(EIGEN3_FOUND)
  remake_find_package(eigen2 CONFIG)
  remake_project_set(EIGEN_VERSION 2 CACHE INTERNAL
    "Version of the eigen linear algebra library")
endif(EIGEN3_FOUND)
remake_find_package(OpenMP OPTIONAL)
remake_find_package(PythonLibs)

if (OPENMP_FOUND)
  remake_define(HAVE_OPENMP ON)
  remake_list_push(CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
endif(OPENMP_FOUND)

remake_include(
  lib
  ${EIGEN${LIBNABO_EIGEN_VERSION}_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIR}
)

remake_add_directories(lib bin)
