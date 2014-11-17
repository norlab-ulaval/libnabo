remake_find_library(ann ANN/ANN.h)
remake_find_package(flann CONFIG)

remake_include(
  ${ANN_HEADERS}/ANN
  ${FLANN_INCLUDE_DIRS}
)
remake_define(HAVE_ANN ON)
remake_define(HAVE_FLANN ON)
  
remake_add_executables(
  knnbench.cpp
  knnbucketsize.cpp
  knnepsilon.cpp
  knnvalidate.cpp
  LINK nabo ${ANN_LIBRARY} flann_cpp
)

remake_test(
  validation validation-2D-exhaustive
  knnvalidate ${CMAKE_SOURCE_DIR}/data/test/scan.2d.txt 10 -1
)
remake_test(
  validation validation-2D-random
  knnvalidate ${CMAKE_SOURCE_DIR}/data/test/scan.2d.txt 10 10000
)
remake_test(
  validation validation-3D-exhaustive
  knnvalidate ${CMAKE_SOURCE_DIR}/data/test/scan.3d.txt 10 -1
)
remake_test(
  validation validation-3D-random
  knnvalidate ${CMAKE_SOURCE_DIR}/data/test/scan.3d.txt 10 10000
)
remake_test(
  validation validation-3D-large-random
  knnvalidate ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 10 1000
)
remake_test(
  validation validation-3D-large-random-radius
  knnvalidate ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 10 1000 0.5
)

remake_test(
  benchmark bench-3D-large-exhaustive-10000-K1
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 1 -10000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-1000-K1
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 1 -1000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-100-K1
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 1 -100 3 5
)
remake_test(
  benchmark bench-3D-large-random-K1
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 1 40000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-10000-K10
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 10 -10000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-1000-K10
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 10 -1000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-100-K10
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 10 -100 3 5
)
remake_test(
  benchmark bench-3D-large-random-K20
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 10 40000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-10000-K30
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 30 -10000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-1000-K30
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 30 -1000 3 5
)
remake_test(
  benchmark bench-3D-large-exhaustive-100-K30
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 30 -100 3 5
)
remake_test(
  benchmark bench-3D-large-random-K30
  knnbench ${CMAKE_SOURCE_DIR}/data/test/scan.3d.large.txt 30 40000 3 5
)
