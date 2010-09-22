include_directories(..)

add_executable(knnvalidate knnvalidate.cpp)
target_link_libraries(knnvalidate nabo)

add_test(validation-2D ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.2d.txt 10 1000)
add_test(validation-3D ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.txt 10 1000)

# find_package(Boost)
# if (Boost_FOUND)
# 	find_path(ANN_INCLUDE_DIR ANN.h
# 		/usr/local/include/ANN
# 		/usr/include/ANN
# 	)
# 	
# 	include_directories(.. ${Boost_INCLUDE_DIRS} ${ANN_INCLUDE_DIR})
# 
# 	add_executable(knnbench knnbench.cpp
# 	# hiiiikkk^2 this is even uglier, let's talk about this irl ok?
# 	~/test/ann_1.1.2/src/ANN.cpp
# 	~/test/ann_1.1.2/src/bd_fix_rad_search.cpp
# 	~/test/ann_1.1.2/src/bd_pr_search.cpp
# 	~/test/ann_1.1.2/src/bd_search.cpp
# 	~/test/ann_1.1.2/src/bd_tree.cpp
# 	~/test/ann_1.1.2/src/brute.cpp
# 	~/test/ann_1.1.2/src/kd_dump.cpp
# 	~/test/ann_1.1.2/src/kd_fix_rad_search.cpp
# 	~/test/ann_1.1.2/src/kd_pr_search.cpp
# 	~/test/ann_1.1.2/src/kd_search.cpp
# 	~/test/ann_1.1.2/src/kd_split.cpp
# 	~/test/ann_1.1.2/src/kd_tree.cpp
# 	~/test/ann_1.1.2/src/kd_util.cpp
# 	~/test/ann_1.1.2/src/perf.cpp)
# 	
# 	target_link_libraries(knnbench nabo)
# else (Boost_FOUND)
# 	message("Disabling benchmark because boost has not been found")
# endif (Boost_FOUND)