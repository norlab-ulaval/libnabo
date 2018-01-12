cmake_minimum_required(VERSION 2.6)

if (NOT CMAKE_VERSION VERSION_LESS "3.1")
	cmake_policy(SET CMP0054 NEW)
endif ()

set(LIB_NAME nabo)
project("lib${LIB_NAME}")

#=============================
# All available CMake options.
#=============================
set(NABO_EIGEN_INCLUDE_DIR "" CACHE STRING "Manually specify the path to Eigen (2 or 3). Empty string means search for Eigen *3*.")
mark_as_advanced(NABO_EIGEN_INCLUDE_DIR)
option(NABO_EXPORT_PACKAGE "Use `export(PACKAGE libnabo)` with CMake?  Makes build directory available for `find_package`, MAY CONFLICT WITH INSTALL." OFF)
mark_as_advanced(NABO_EXPORT_PACKAGE)
set(NABO_DOC_INSTALL_TARGET "share/doc/${PROJECT_NAME}/api" CACHE STRING "Target where to install doxygen documentation.")
option(NABO_SHARED_LIBS "Set to ON to build shared library (default: OFF)." OFF)
# OpenCL disabled as its code is not up-to-date with API
option(NABO_USE_OPEN_CL "Set to ON to look for OpenCL (note: OpenCL code not up to date with API)" OFF)
option(NABO_USE_OPEN_MP "Use OpenMP by default (if OpenMP is available)? Set to OFF to not use OpenMP." ON)

# end options

# Extract version from header
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
file(
	STRINGS                       # Read in a file to memory.
	nabo/nabo.h                   # File to parse for version number.
	NABO_PROJECT_VERSION          # Where to store the results (should just be one line)
	REGEX "#define NABO_VERSION " # The space after is important to ignore 'NABO_VERSION_INT'
)
# If no matches found, something is wrong with nabo.h
if (NOT NABO_PROJECT_VERSION)
	message(SEND_ERROR "Cannot find version number in '${CMAKE_CURRENT_SOURCE_DIR}/nabo/nabo.h'.")
endif (NOT NABO_PROJECT_VERSION)
# Transform '#define NABO_VERSION "X.Y.Z"' into 'X.Y.Z'
string(REGEX REPLACE ".*\"(.*)\".*" "\\1" NABO_PROJECT_VERSION "${NABO_PROJECT_VERSION}" )

if (NOT CMAKE_BUILD_TYPE)
	message("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
	set(CMAKE_BUILD_TYPE Release CACHE STRING
	  "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	  FORCE)
else (NOT CMAKE_BUILD_TYPE)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("\n=================================================================================")
		message("\n-- Build type: Debug. Performance will be terrible!")
		message("-- Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build.")
		message("\n=================================================================================")
	endif (CMAKE_BUILD_TYPE STREQUAL "Debug")
endif (NOT CMAKE_BUILD_TYPE)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-O3)
endif(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")

# Documentation
set(DOXYFILE_LATEX false)
include(UseDoxygen)

# Switch on warnings.
add_definitions(-Wall -Wextra)

# Check the compiler version as we need C++11 support.
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	# using Clang
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
		message(FATAL_ERROR "Your clang compiler has version ${CMAKE_CXX_COMPILER_VERSION}, while version 3.3 or later is required")
	endif ()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	# using AppleClang
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.1")
		message(FATAL_ERROR "Your XCode environment has version ${CMAKE_CXX_COMPILER_VERSION}, while version 5.1 or later is required")
	endif ()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	# using GCC
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.2")
		message(FATAL_ERROR "Your g++ compiler has version ${CMAKE_CXX_COMPILER_VERSION}, while version 4.8.2 or later is required")
	endif ()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# using MSVC
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.0.23506")
		message(FATAL_ERROR "Your Microsoft Visual C++ compiler has version ${CMAKE_CXX_COMPILER_VERSION}, while version MSVC 2015 Update 1+ is required")
	endif ()
endif ()

# enable C++11 support.
if (CMAKE_VERSION VERSION_LESS "3.1")
	if (MSVC)
		message(FATAL_ERROR "CMake version 3.1 or later is required to compile ${PROJECT_NAME} with Microsoft Visual C++")
	endif ()
	if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		set (CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
	else ()
		set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
	endif ()
else ()
	set (CMAKE_CXX_STANDARD 11)
endif ()

# Boost required for any and cstdint on MSVC < 2010
find_package(Boost COMPONENTS REQUIRED)
if(MSVC AND (MSVC_VERSION LESS 1600))
	add_definitions(-DBOOST_STDINT)
endif(MSVC AND (MSVC_VERSION LESS 1600))

# openmp
if (NABO_USE_OPEN_MP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		message(STATUS "Nabo: OpenMP found!")
		add_definitions(-fopenmp -DHAVE_OPENMP)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		if (CMAKE_COMPILER_IS_GNUCC)
			set(EXTRA_LIBS ${EXTRA_LIBS} gomp)
		endif(CMAKE_COMPILER_IS_GNUCC)
	else (OPENMP_FOUND)
		message(STATUS "Nabo: OpenMP *NOT* found -- skipping.")
	endif(OPENMP_FOUND)
endif (NABO_USE_OPEN_MP)

# eigen 2 or 3 supported, if 'NABO_EIGEN_INCLUDE_DIR' not specified then only v3 is
# searched for.
if (NABO_EIGEN_INCLUDE_DIR)
	# If specified, make sure Eigen/Core exists.
	if (NOT EXISTS "${NABO_EIGEN_INCLUDE_DIR}/Eigen/Core")
		message(FATAL_ERROR "The file '${NABO_EIGEN_INCLUDE_DIR}/Eigen/Core' does not exist!")
	endif (NOT EXISTS "${NABO_EIGEN_INCLUDE_DIR}/Eigen/Core")
else()
	if (${CMAKE_VERSION} VERSION_LESS 3.0.0)
		find_package(PkgConfig)
		if (NOT PkgConfig_FOUND)
			message(FATAL_ERROR "PkgConfig not available, please specify 'NABO_EIGEN_INCLUDE_DIR' manually.")
		endif (NOT PkgConfig_FOUND)
		pkg_search_module(Eigen3 REQUIRED eigen3)
		set(NABO_EIGEN_INCLUDE_DIR "${Eigen3_INCLUDE_DIRS}")
	else (${CMAKE_VERSION} VERSION_LESS 3.0.0)
		find_package(Eigen3 REQUIRED)# FindEigen3.cmake introduced in 3.0
		set(NABO_EIGEN_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR}")
	endif (${CMAKE_VERSION} VERSION_LESS 3.0.0)
endif()


# optionally, opencl
if (NABO_USE_OPEN_CL)
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
		message("OpenCL enabled and found, enabling CL support")
	else (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)
		message("OpenCL enabled but not found, disabling CL support")
	endif (OPENCL_INCLUDE_DIR AND OPENCL_LIBRARIES)
else (NABO_USE_OPEN_CL)
	message("OpenCL disabled, not looking for it")
endif (NABO_USE_OPEN_CL)

# include all libs so far
include_directories(${CMAKE_SOURCE_DIR} ${NABO_EIGEN_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

# main nabo lib
set(NABO_SRC
	nabo/nabo.cpp
	nabo/brute_force_cpu.cpp
	nabo/kdtree_cpu.cpp
	nabo/kdtree_opencl.cpp
)

if (NABO_SHARED_LIBS)
	add_library(${LIB_NAME} SHARED ${NABO_SRC})
	install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib)
else (NABO_SHARED_LIBS)
	add_library(${LIB_NAME} ${NABO_SRC})
	add_definitions(-fPIC)
	install(TARGETS ${LIB_NAME} ARCHIVE DESTINATION lib)
endif (NABO_SHARED_LIBS)
set_target_properties(${LIB_NAME} PROPERTIES VERSION "${NABO_PROJECT_VERSION}" SOVERSION 1)

# create doc before installing
install(FILES nabo/nabo.h DESTINATION include/nabo)
install(FILES nabo/third_party/any.hpp DESTINATION include/nabo/third_party)
install(FILES README.md DESTINATION share/doc/${PROJECT_NAME})
if (DOXYGEN_FOUND)
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION ${NABO_DOC_INSTALL_TARGET})
endif (DOXYGEN_FOUND)

enable_testing()

add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(python)

# Install catkin package.xml
install(FILES package.xml DESTINATION share/libnabo)

#=============================================
# to allow find_package() on libnabo
#=============================================
#
# the following case be used in an external project requiring libnabo:
#  ...
#  find_package(libnabo)
#  include_directories(${libnabo_INCLUDE_DIRS})
#  target_link_libraries(executableName ${libnabo_LIBRARIES})
#  ...

# NOTE: the following will support find_package for 1) local build (make) and 2) for installed files (make install)

# 1- local build #

# Register the local build in case one doesn't use "make install"
if (NABO_EXPORT_PACKAGE)
	export(PACKAGE libnabo)
endif (NABO_EXPORT_PACKAGE)

# Create variable with the library location
set(libnabo_library $<TARGET_FILE:${LIB_NAME}>)

# Create variable for the local build tree
get_property(libnabo_include_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# Configure & generate config file for local build tree
configure_file(libnaboConfig.cmake.in
	"${PROJECT_BINARY_DIR}/libnaboConfig.cmake.conf" @ONLY)
file(GENERATE
	OUTPUT "${PROJECT_BINARY_DIR}/libnaboConfig.cmake"
	INPUT "${PROJECT_BINARY_DIR}/libnaboConfig.cmake.conf")

# 2- installation build #

# Change the library location for an install location
set(libnabo_library ${CMAKE_INSTALL_PREFIX}/lib/$<TARGET_FILE_NAME:${LIB_NAME}>)

# Change the include location for the case of an install location
set(libnabo_include_dirs ${CMAKE_INSTALL_PREFIX}/include)

# We put the generated file for installation in a different repository (i.e., ./CMakeFiles/)
configure_file(libnaboConfig.cmake.in
	"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libnaboConfig.cmake.conf" @ONLY)
file(GENERATE
	OUTPUT "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libnaboConfig.cmake"
	INPUT "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libnaboConfig.cmake.conf")

# The same versioning file can be used for both cases
configure_file(libnaboConfigVersion.cmake.in
	"${PROJECT_BINARY_DIR}/libnaboConfigVersion.cmake" @ONLY)

install(FILES
	"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libnaboConfig.cmake"
	"${PROJECT_BINARY_DIR}/libnaboConfigVersion.cmake"
	DESTINATION share/libnabo/cmake COMPONENT dev)


#=============================================
# Add uninstall target
#=============================================
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

add_custom_target(uninstall
	COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
