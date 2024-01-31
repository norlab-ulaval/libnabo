#!/bin/bash
# =================================================================================================
#
# Install all libnabo dependencies
#
# Redirect the execution to libnabo-build-system library
#
# Usage:
#   $ bash libnabo_dependencies_installer.bash [--test-run]
#
# =================================================================================================
PARAMS="$@"

MSG_DIMMED_FORMAT="\033[1;2m"
MSG_ERROR_FORMAT="\033[1;31m"
MSG_END_FORMAT="\033[0m"

function nabo::install_libnabo_dependencise(){

  # ....path resolution logic......................................................................
  NABO_ROOT="$(dirname "$(realpath "$0")")"

  cd "${NABO_ROOT}" || exit 1

  # ....Load environment variables from file.......................................................
  # . . Source NABO environment variables  . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  if [[ -f .env.libnabo ]]; then
    set -o allexport && source .env.libnabo && set +o allexport
  else
    echo -e "${MSG_ERROR_FORMAT}[NABO ERROR]${MSG_END_FORMAT} .env.libnabo unreachable. Cwd $(pwd)" 1>&2
    exit 1
  fi

  # . . Source NBS dependencies . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  cd "${N2ST_PATH}" || exit 1
  source import_norlab_shell_script_tools_lib.bash

  # . . Source NABO-build-system environment variables . . . . . . . . . . . . . . . . . . . . . ..
  cd "${NABO_BUILD_SYSTEM_PATH}" || exit 1
  if [[ -f .env ]]; then
    set -o allexport && source .env && set +o allexport
  else
    echo -e "${MSG_ERROR_FORMAT}[NABO ERROR]${MSG_END_FORMAT} .env unreachable. Cwd $(pwd)" 1>&2
    exit 1
  fi


  # ====Begin======================================================================================
  norlab_splash "${NBS_SPLASH_NAME}" "https://github.com/${NBS_REPOSITORY_DOMAIN:?err}/${NBS_REPOSITORY_NAME:?err}"
  export SHOW_SPLASH_IDU=false
  export SHOW_SPLASH_IDDU=false

  # ....Install general dependencies...............................................................
  cd "${NABO_PATH:?err}"/build_system/ubuntu || exit 1
  # shellcheck disable=SC2086
  source nabo_install_dependencies_ubuntu.bash ${PARAMS[@]}

  # ....Install documentation related dependencies.................................................
  cd "${NABO_PATH:?err}"/build_system/ubuntu || exit 1
  # shellcheck disable=SC2086
  source nabo_install_doc_dependencies_ubuntu.bash ${PARAMS[@]}

  print_msg_done "All Libnabo dependencies installed"
}

# ::::Main:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if [[ "${BASH_SOURCE[0]}" = "$0" ]]; then
  # This script is being run, ie: __name__="__main__"
  nabo::install_libnabo_dependencise
else
  # This script is being sourced, ie: __name__="__source__"
  echo -e "${MSG_ERROR_FORMAT}[NABO ERROR]${MSG_END_FORMAT} Execute this script in a subshell i.e.: $ bash libnabo_dependencies_installer.bash" 1>&2
  exit 1
fi
