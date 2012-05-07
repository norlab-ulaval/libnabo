cmake_minimum_required(VERSION 2.6)

set(LIB_NAME nabo)
project("lib${LIB_NAME}")

# Extract version from header
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
execute_process(
	COMMAND grep "NABO_VERSION " nabo/nabo.h
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	OUTPUT_VARIABLE PROJECT_VERSION
	OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REGEX REPLACE ".*\"(.*)\".*" "\\1" PROJECT_VERSION "${PROJECT_VERSION}" )

# Documentation
set(DOXYFILE_LATEX false)
include(UseDoxygen)

# Compilation
add_definitions(-Wall)

# Boost required for any and cstdint on MSVC < 2010
find_package(Boost COMPONENTS REQUIRED)
if(MSVC AND (MSVC_VERSION LESS 1600))
	add_definitions(-DBOOST_STDINT)
endif(MSVC AND (MSVC_VERSION LESS 1600))

# openmp
# find_package(OpenMP)
# if (OPENMP_FOUND)
# 	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
# endif(OPENMP_FOUND)

# eigen 2 or 3
find_path(EIGEN_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen3
	/usr/local/include/eigen2
	/usr/local/include/eigen
	/usr/include/eigen3
	/usr/include/eigen2
	/usr/include/eigen
	/opt/ros/diamondback/stacks/geometry/eigen/include
)

# optionally, opencl
# OpenCL disabled as its code is not up-to-date with API
# set(USE_OPEN_CL "false" CACHE BOOL "Set to ON to look for OpenCL")
# if (USE_OPEN_CL)
# 	find_path(OPENCL_INCLUDE_DIR CL/cl.h
# 		/usr/local/include
# 		/usr/include
# 	)
# 	if (WIN32)
# 		find_library(OPENCL_LIBRARIES opencl64)
# 		if (!OPENCL_LIBRARIES)
# 			find_library(OPENCL_LIBRARIES opencl32)
# 		endif (!OPENCL_LIBRARIES)
# 	else (WIN32)
# 		find_library(OPENCL_LIBRARIES OpenCL ENV LD_LIBRARY_PATH)
# 	endif (WIN32)
# 	if (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)
# 		add_definitions(-DHAVE_OPENCL)
# 		set(EXTRA_LIBS ${OPENCL_LIBRARIES} ${EXTRA_LIBS})
# 		include_directories(${OPENCL_INCLUDE_DIR})
# 		add_definitions(-DOPENCL_SOURCE_DIR=\"${CMAKE_SOURCE_DIR}/nabo/opencl/\")
# 		message("OpenCL enabled and found, enabling CL support")
# 	else (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)
# 		message("OpenCL enabled but not found, disabling CL support")
# 	endif (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)
# else(USE_OPEN_CL)
# 	message("OpenCL disabled, not looking for it")
# endif(USE_OPEN_CL)

# include all libs so far
include_directories(${EIGEN_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

# main nabo lib
set(NABO_SRC
	nabo/nabo.cpp
	nabo/brute_force_cpu.cpp
	nabo/kdtree_cpu.cpp
	nabo/kdtree_opencl.cpp
)
set(SHARED_LIBS "false" CACHE BOOL "To build shared (true) or static (false) library")
if (SHARED_LIBS)
	add_library(${LIB_NAME} SHARED ${NABO_SRC})
	install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib)
else(SHARED_LIBS)
	add_library(${LIB_NAME} ${NABO_SRC})
	install(TARGETS ${LIB_NAME} ARCHIVE DESTINATION lib)
endif(SHARED_LIBS)
set_target_properties(${LIB_NAME} PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1)

# create doc before installing
set(DOC_INSTALL_TARGET "share/doc/${PROJECT_NAME}/api" CACHE STRING "Target where to install doxygen documentation")
add_dependencies(${LIB_NAME} doc)
install(FILES nabo/nabo.h DESTINATION include/nabo)
install(FILES README.md DESTINATION share/doc/${PROJECT_NAME})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION ${DOC_INSTALL_TARGET})

enable_testing()

add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(python)