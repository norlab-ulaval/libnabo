#!/bin/bash
#
# Utility script to create a multi-architecture docker builder
#
# usage:
#   $ bash ./nabo_utility_script/nabo_create_multiarch_docker_builder.bash
#
function nabo::create_multiarch_docker_builder() {
  local TMP_CWD
  TMP_CWD=$(pwd)

  # ....Project root logic.........................................................................
  NABO_PATH=$(git rev-parse --show-toplevel)

  # ====Begin=====================================================================================
  cd "${NABO_PATH}"/build_system/utilities/norlab-build-system/install_scripts/ \
    && bash nbs_create_multiarch_docker_builder.bash

  # ====Teardown===================================================================================
  cd "${TMP_CWD}"
}

# ::::main:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
nabo::create_multiarch_docker_builder
