cmake_minimum_required(VERSION 2.6)

project(nabo)

# eigen
find_path(EIGEN2_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen2
	/usr/include/eigen2
)

# include all libs so far
include_directories(${EIGEN2_INCLUDE_DIR})

# our lib
add_library(nabo nabo.cpp brute_force_cpu.cpp kdtree_cpu.cpp)

# our test program
add_executable(knn knn.cpp)
target_link_libraries(knn nabo)