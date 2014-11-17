cmake_minimum_required(VERSION 2.6)

include(ReMake)

remake_project(
  libnabo
  VERSION 1.0.4
  RELEASE stable
  SUMMARY "Fast k-Nearest Neighbor library for low-dimensional spaces."
  AUTHOR "Stephane Magnenat"
  CONTACT "stephane.magnenat@mavt.ethz.ch"
  HOME "http://github.com/ethz-asl/libnabo"
  LICENSE "BSD 3-Clause License (BSD-3-Clause)"
  README README.md
  INSTALL /usr
  HEADER_DESTINATION include/nabo
  PREFIX nabo-
)

remake_add_directories(
  python
  data
)
remake_doc(
  html
)

remake_pack_deb(
  DEPENDS libgomp1
)
remake_pack_deb(
  COMPONENT python
  DEPENDS libnabo python-numpy libpython2.[0-9]* libboost-python[0-9.]*
  DESCRIPTION "Python bindings"
)
remake_pack_deb(
  COMPONENT test
  DEPENDS libnabo libnabo-python libann0 libflann[0-9.]*
  DESCRIPTION "testing"
)
remake_pack_deb(
  COMPONENT dev
  DEPENDS libnabo libeigen${LIBNABO_EIGEN_VERSION}-dev
  DESCRIPTION "development headers"
)
remake_pack_deb(
  ARCH all
  COMPONENT doc
  DESCRIPTION "documentation"
)

remake_distribute_deb(
  DISTRIBUTION lucid
  SECTION libs
  UPLOAD ppa:kralf/asl
  DEPENDS python-dev libgomp1 libboost-python-dev libeigen3-dev python-numpy
    libann-dev libflann-dev remake pkg-config doxygen
  PASS CMAKE_BUILD_TYPE LIBNABO_GIT_REVISION
)
remake_distribute_deb(
  DISTRIBUTION precise
  SECTION libs
  UPLOAD ppa:kralf/asl
  DEPENDS python-dev libgomp1 libboost-python-dev libeigen3-dev python-numpy
    libann-dev libflann-dev remake pkg-config doxygen
  PASS CMAKE_BUILD_TYPE LIBNABO_GIT_REVISION
)
remake_distribute_deb(
  DISTRIBUTION trusty
  SECTION libs
  UPLOAD ppa:kralf/asl
  DEPENDS python-dev libgomp1 libboost-python-dev libeigen3-dev python-numpy
    libann-dev libflann-dev remake pkg-config doxygen
  PASS CMAKE_BUILD_TYPE LIBNABO_GIT_REVISION
)
