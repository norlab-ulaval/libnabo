#!/bin/bash
#
# Docker entrypoint for development task
#
# Usage:
#   $ bash entrypoint.bash [<any-cmd>]
#
# Parameter
#   <any-cmd>      Optional command executed in a subprocess at the end of the entrypoint script.
#

# ....Load environment variables from file.........................................................
set -o allexport
source "${NBS_LIB_INSTALL_PATH}/${NBS_REPOSITORY_NAME}/build_system/.env"
set +o allexport

# ....Helper function..............................................................................
source "${N2ST_PATH:?err}/import_norlab_shell_script_tools_lib.bash"

# ====Begin========================================================================================
clear

norlab_splash "${NBS_SPLASH_NAME:?err}" "https://github.com/${NBS_REPOSITORY_DOMAIN:?err}/${NBS_REPOSITORY_NAME:?err}"

#ifconfig
pwd && tree -L 1

# ====Continue=====================================================================================================
exec "$@"
