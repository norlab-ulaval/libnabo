#!/bin/bash

# ....path resolution logic........................................................................
_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[0]}")"
NABO_ROOT_DIR="$(dirname "${_PATH_TO_SCRIPT}")"
cd "${NABO_ROOT_DIR}/../../"

# ====begin========================================================================================
#bash nabo_crawl_libnabo_build_matrix.bash --fail-fast -- build --dry-run ci_PR_amd64 ci_PR_arm64v8
bash nabo_crawl_libnabo_build_matrix.bash --fail-fast -- build --dry-run
