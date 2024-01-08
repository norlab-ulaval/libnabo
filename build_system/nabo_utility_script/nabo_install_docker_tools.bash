#!/bin/bash
#
# Utility script to install docker related tools and execute basic configuration
#
# usage:
#   $ bash ./nabo_utility_script/nabo_install_docker_tools.bash
#

function nabo::install_docker_tools() {
  local TMP_CWD
  TMP_CWD=$(pwd)

  # ....Project root logic.........................................................................
  NABO_PATH=$(git rev-parse --show-toplevel)

  # ====Begin=====================================================================================
  cd "${NABO_PATH}"/build_system/utilities/norlab-build-system/install_scripts/ \
    && bash nbs_install_docker_tools.bash

  # ====Teardown===================================================================================
  cd "${TMP_CWD}"
}

# ::::main:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
nabo::install_docker_tools
