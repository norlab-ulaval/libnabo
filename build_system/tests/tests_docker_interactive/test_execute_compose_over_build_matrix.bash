#!/bin/bash

clear

# ....path resolution logic........................................................................
_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[0]}")"
NABO_ROOT_DIR="$(dirname "${_PATH_TO_SCRIPT}")"
cd "${NABO_ROOT_DIR}/../../.."

pwd

# ====begin========================================================================================
bash build_system/tests/tests_docker_interactive/build_and_run_IamBuildSystemTester.bash "bash ./nabo_crawl_libnabo_build_matrix.bash \
    --repository-version-build-matrix-override latest  \
    --cmake-build-type-build-matrix-override None \
    --docker-debug-logs \
    -- build --dry-run"
