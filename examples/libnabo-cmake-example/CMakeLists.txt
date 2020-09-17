cmake_minimum_required(VERSION 3.0)
project(libnabo-cmake-example)

# This looks for the libnabo-config.cmake file, which in turns
# includes libnabo-targets.cmake:
find_package(libnabo REQUIRED)

add_executable(${PROJECT_NAME} trivial.cpp)

# This sets both the include directories and link library:
target_link_libraries(${PROJECT_NAME} libnabo::nabo)
