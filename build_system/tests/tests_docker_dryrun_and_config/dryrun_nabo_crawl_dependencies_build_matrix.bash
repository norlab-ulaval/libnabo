#!/bin/bash

# ....path resolution logic........................................................................
_PATH_TO_SCRIPT="$(realpath "${BASH_SOURCE[0]}")"
NABO_ROOT_DIR="$(dirname "${_PATH_TO_SCRIPT}")"
cd "${NABO_ROOT_DIR}/../../"

# ====begin========================================================================================
#FLAGS=( '--fail-fast' '--' 'build' '--dry-run' )
#bash nabo_crawl_dependencies_build_matrix.bash ${FLAGS[@]}

#bash nabo_crawl_dependencies_build_matrix.bash --fail-fast -- build --dry-run dependencies dependencies-doc

bash nabo_crawl_dependencies_build_matrix.bash --fail-fast -- build --dry-run
