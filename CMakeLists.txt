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
add_library(nabo 
	nabo/nabo.cpp
	nabo/brute_force_cpu.cpp
	nabo/kdtree_cpu.cpp
)
install(TARGETS nabo ARCHIVE DESTINATION lib)
install(FILES nabo/nabo.h DESTINATION include/nabo)

add_subdirectory(benchmarks)
add_subdirectory(tests)

