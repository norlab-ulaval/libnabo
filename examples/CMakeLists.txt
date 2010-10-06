include_directories(..)

add_executable(trivial trivial.cpp)
target_link_libraries(trivial nabo ${Boost_LIBRARIES})