include (CheckSymbolExists) 

include_directories(..)

add_executable(knnvalidate knnvalidate.cpp)
target_link_libraries(knnvalidate nabo  ${EXTRA_LIBS} ${Boost_LIBRARIES})
# disable as experimental is broken target_link_libraries(knnvalidate naboexperimental nabo  ${EXTRA_LIBS} ${Boost_LIBRARIES})

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
	set(EXTRA_LIBS ${ANN_LIBRARY} ${EXTRA_LIBS})
	add_definitions(-DHAVE_ANN)
else (ANN_INCLUDE_DIR AND ANN_LIBRARY)
	message("ANN not found, disabling it in benchmarks")
endif (ANN_INCLUDE_DIR AND ANN_LIBRARY)

find_path(FLANN_INCLUDE_DIR flann/flann.hpp
	/usr/local/include/
	/usr/include
)
find_library(FLANN_LIBRARY
	flann_cpp
	PATHS /usr/local/lib
	/usr/lib
)

if (FLANN_INCLUDE_DIR AND FLANN_LIBRARY)
	include_directories(${FLANN_INCLUDE_DIR})
	set(EXTRA_LIBS ${FLANN_LIBRARY} ${EXTRA_LIBS})
	add_definitions(-DHAVE_FLANN)
else (FLANN_INCLUDE_DIR AND FLANN_LIBRARY)
	message("FLANN not found, disabling it in benchmarks")
endif (FLANN_INCLUDE_DIR AND FLANN_LIBRARY)

check_symbol_exists(_POSIX_TIMERS "unistd.h;time.h" POSIX_TIMERS)
if (POSIX_TIMERS)
	set(EXTRA_LIBS ${EXTRA_LIBS} rt)
endif (POSIX_TIMERS)

include_directories(..)

add_executable(knnbench knnbench.cpp)
target_link_libraries(knnbench nabo ${EXTRA_LIBS} ${Boost_LIBRARIES})
# disable as experimental is broken target_link_libraries(knnbench naboexperimental nabo ${EXTRA_LIBS} ${Boost_LIBRARIES})

add_test(bench-3D-large-exhaustive-10000-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -10000 3 5)
add_test(bench-3D-large-exhaustive-1000-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -1000 3 5)
add_test(bench-3D-large-exhaustive-100-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 -100 3 5)
add_test(bench-3D-large-random-K10 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 10 40000 3 5)
add_test(bench-3D-large-exhaustive-10000-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 -10000 3 5)
add_test(bench-3D-large-exhaustive-1000-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 -1000 3 5)
add_test(bench-3D-large-exhaustive-100-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 -100 3 5)
add_test(bench-3D-large-random-K20 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 20 40000 3 5)
add_test(bench-3D-large-exhaustive-10000-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 30 -10000 3 5)
add_test(bench-3D-large-exhaustive-1000-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 30 -1000 3 5)
add_test(bench-3D-large-exhaustive-100-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 30 -100 3 5)
add_test(bench-3D-large-random-K100 ${EXECUTABLE_OUTPUT_PATH}/knnbench ${CMAKE_CURRENT_SOURCE_DIR}/data/scan.3d.large.txt 30 40000 3 5)

add_executable(knnepsilon knnepsilon.cpp)
target_link_libraries(knnepsilon nabo ${EXTRA_LIBS} ${Boost_LIBRARIES})