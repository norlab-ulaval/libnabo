remake_add_library(nabo PREFIX OFF)
remake_add_headers(nabo.h)
remake_pkg_config_generate(REQUIRES eigen${LIBNABO_EIGEN_VERSION})
