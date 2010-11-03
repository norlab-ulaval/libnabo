cmake_minimum_required(VERSION 2.6)

project(libnabo)

# Documentation, extract version from header
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
execute_process(
	COMMAND grep "NABO_VERSION " nabo/nabo.h
	COMMAND sed "s/.*\"\(.*\)\"/\1/g"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	OUTPUT_VARIABLE PROJECT_VERSION
	OUTPUT_STRIP_TRAILING_WHITESPACE
)
include(UseDoxygen)

# Compilation
add_definitions(-Wall)

find_package(Boost REQUIRED)

# C++0x 
if(CMAKE_COMPILER_IS_GNUCXX) #Maybe "OR MINGW"
	add_definitions(-std=c++0x) #Enable C++0x mode
endif(CMAKE_COMPILER_IS_GNUCXX)
# openmp
# find_package(OpenMP)
# if (OPENMP_FOUND)
# 	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
# endif(OPENMP_FOUND)

# eigen 2
find_path(EIGEN2_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen2
	/usr/include/eigen2
)

# opencl
find_path(OPENCL_INCLUDE_DIR CL/cl.h
	/usr/local/include
	/usr/include
)
if (WIN32)
	find_library(OPENCL_LIBRARIES opencl64)
	if (!OPENCL_LIBRARIES)
		find_library(OPENCL_LIBRARIES opencl32)
	endif (!OPENCL_LIBRARIES)
else (WIN32)
	find_library(OPENCL_LIBRARIES OpenCL ENV LD_LIBRARY_PATH)
endif (WIN32)
if (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)
	add_definitions(-DHAVE_OPENCL)
	set(EXTRA_LIBS ${OPENCL_LIBRARIES} ${EXTRA_LIBS})
	include_directories(${OPENCL_INCLUDE_DIR})
  add_definitions(-DOPENCL_SOURCE_DIR=\"${CMAKE_SOURCE_DIR}/nabo/opencl/\")
endif (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)

# include all libs so far
include_directories(${EIGEN2_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

# main nabo lib
add_library(nabo 
	nabo/nabo.cpp
	nabo/brute_force_cpu.cpp
	nabo/kdtree_cpu.cpp
	nabo/kdtree_opencl.cpp
)
install(TARGETS nabo ARCHIVE DESTINATION lib)
install(FILES nabo/nabo.h DESTINATION include/nabo)

enable_testing()
add_subdirectory(examples)
#add_subdirectory(experimental)
add_subdirectory(tests)

