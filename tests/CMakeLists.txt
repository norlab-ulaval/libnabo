include_directories(..)

add_executable(knnvalidate knnvalidate.cpp)
target_link_libraries(knnvalidate nabo)

add_test(validation-2D-exhaustive ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.2d.txt 10 -1)
add_test(validation-2D-random ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.2d.txt 10 10000)
add_test(validation-3D-exhaustive ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.txt 10 -1)
add_test(validation-3D-random ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.txt 10 10000)
add_test(validation-3D-large-random ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 1000)

find_package(Boost)
find_path(ANN_INCLUDE_DIR ANN.h
	/usr/local/include/ANN
	/usr/include/ANN
)
find_library(ANN_LIBRARY
	NAMES ANN ann
	PATHS /usr/local/lib
	/usr/lib
)

if (Boost_FOUND AND ANN_INCLUDE_DIR AND ANN_LIBRARY)
	include_directories(.. ${Boost_INCLUDE_DIRS} ${ANN_INCLUDE_DIR})

	add_executable(knnbench knnbench.cpp)
	target_link_libraries(knnbench nabo ${ANN_LIBRARY})
	
	add_test(bench-3D-large-exhaustive ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -1)
	add_test(bench-3D-large-random ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 10000)
else (Boost_FOUND AND ANN_INCLUDE_DIR AND ANN_LIBRARY)
	message("Disabling performance benchmark because boost or Ann has not been found")
endif (Boost_FOUND AND ANN_INCLUDE_DIR AND ANN_LIBRARY)