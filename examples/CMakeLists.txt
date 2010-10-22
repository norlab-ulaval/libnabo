include_directories(..)

add_executable(trivial trivial.cpp)
target_link_libraries(trivial nabo ${EXTRA_LIBS} ${Boost_LIBRARIES})

add_executable(usage usage.cpp)
target_link_libraries(usage nabo ${EXTRA_LIBS} ${Boost_LIBRARIES})