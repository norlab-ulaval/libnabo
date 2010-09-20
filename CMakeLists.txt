project(knn)

# eigen
find_path(EIGEN2_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen2
	/usr/include/eigen2
)

# include all libs so far
include_directories(${EIGEN2_INCLUDE_DIR} )

# our lib
add_library(nabo knn.cpp)

# our test program
add_executable(knn knn.cpp)