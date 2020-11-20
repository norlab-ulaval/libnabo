# Set compiler flags
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-O3)
endif(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")

# Find catkin macros and libraries
find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS python)

find_package(OpenMP)
if (NOT OpenMP_FOUND)
  message("OpenMP was not found. It is highly recommended to build libnabo with OpenMP support.")
endif()

# Catkin package macro
catkin_package(
  INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
  LIBRARIES
    nabo
)

########################
## Library definition ##
########################
# Nabo
add_library(nabo
  ${NABO_SRC}
  ${NABO_HEADERS}
)
target_include_directories(nabo PUBLIC
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/nabo
)
target_include_directories(nabo SYSTEM PRIVATE
  ${EIGEN3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${OpenMP_CXX_INCLUDE_DIRS}
)
target_link_libraries(nabo PRIVATE
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  ${OpenMP_CXX_LIBRARIES}
)
target_compile_options(nabo PRIVATE
  "${OpenMP_CXX_FLAGS}"
)
target_compile_definitions(nabo PRIVATE
  -DHAVE_OPENMP=${OpenMP_FOUND}
)

# Python bindings
add_subdirectory(python)

#############
## Install ##
#############
install(
  TARGETS
    nabo
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(
  DIRECTORY
    ${CMAKE_SOURCE_DIR}/nabo/
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}/nabo
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
)