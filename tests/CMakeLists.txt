include_directories(..)

add_executable(knnvalidate knnvalidate.cpp)
target_link_libraries(knnvalidate nabo ${Boost_LIBRARIES})


add_test(validation-2D-exhaustive ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.2d.txt 10 -1)
add_test(validation-2D-random ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.2d.txt 10 10000)
add_test(validation-3D-exhaustive ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.txt 10 -1)
add_test(validation-3D-random ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.txt 10 10000)
add_test(validation-3D-large-random ${EXECUTABLE_OUTPUT_PATH}/knnvalidate ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 1000)

find_path(ANN_INCLUDE_DIR ANN.h
	/usr/local/include/ANN
	/usr/include/ANN
)
find_library(ANN_LIBRARY
	NAMES ANN ann
	PATHS /usr/local/lib
	/usr/lib
)

if (ANN_INCLUDE_DIR AND ANN_LIBRARY)
	include_directories(${ANN_INCLUDE_DIR})
	set(EXTRA_LIBS ${ANN_LIBRARY})
	add_definitions(-DHAVE_ANN)
else (ANN_INCLUDE_DIR AND ANN_LIBRARY)
	message("ANN not found, disabling it in benchmarks")
endif (ANN_INCLUDE_DIR AND ANN_LIBRARY)

include_directories(.. ${Boost_INCLUDE_DIRS})

add_executable(knnbench knnbench.cpp)
target_link_libraries(knnbench nabo ${EXTRA_LIBS} ${Boost_LIBRARIES})

add_test(bench-3D-large-exhaustive-10000-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -10000 3)
add_test(bench-3D-large-exhaustive-1000-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -1000 3)
add_test(bench-3D-large-exhaustive-100-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -100 3)
add_test(bench-3D-large-random-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 40000 3)
add_test(bench-3D-large-exhaustive-10000-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 -10000 3)
add_test(bench-3D-large-exhaustive-1000-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 -1000 3)
add_test(bench-3D-large-exhaustive-100-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 -100 3)
add_test(bench-3D-large-random-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 40000 3)
add_test(bench-3D-large-exhaustive-10000-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 100 -10000 3)
add_test(bench-3D-large-exhaustive-1000-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 100 -1000 3)
add_test(bench-3D-large-exhaustive-100-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 100 -100 3)
add_test(bench-3D-large-random-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 100 40000 3)
