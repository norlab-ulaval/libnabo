cmake_minimum_required(VERSION 3.5.1)
project(libnabo)

# Source and header files.
set(NABO_SRC
  nabo/nabo.cpp
  nabo/brute_force_cpu.cpp
  nabo/kdtree_cpu.cpp
  nabo/kdtree_opencl.cpp
)

file(GLOB_RECURSE NABO_HEADERS "nabo/*.h")

# Global compiler flags.
add_compile_options(-Wall -Wextra -Wpedantic)

# Build type options:
#	* catkin.
#	* cmake.
#	* ament (future).
set(BUILD_TYPE "catkin")
if (BUILD_TYPE STREQUAL "catkin")
  include(cmake/CatkinBuild.cmake)
else()
  # For a pure-cmake build, follow the instructions on the package.xml file, which
  # specify how to set the buildtool to cmake and export the package accordingly.
  include(cmake/PureCMakeBuild.cmake)
endif()
