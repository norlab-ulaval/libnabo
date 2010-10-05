cmake_minimum_required(VERSION 2.6)

project(nabo)

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

# include all libs so far
include_directories(${EIGEN2_INCLUDE_DIR})

# our lib
add_library(nabo 
	nabo/nabo.cpp
	nabo/brute_force_cpu.cpp
	nabo/kdtree_cpu.cpp
)
install(TARGETS nabo ARCHIVE DESTINATION lib)
install(FILES nabo/nabo.h DESTINATION include/nabo)

enable_testing()
add_subdirectory(tests)

