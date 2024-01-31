#!/bin/bash
# =================================================================================================
#
# Run Libnabo unit test
#
# Usage:
#   $ bash nabo_execute_nabo_unittest.bash
#
# Notes:
#   The script propagate the utest exit code on exit
# =================================================================================================

# ....Project root logic...........................................................................
TMP_CWD=$(pwd)

NABO_PATH=$(git rev-parse --show-toplevel)
cd "${NABO_PATH}/build_system" || exit 1

# ....Load environment variables from file.........................................................
set -o allexport
source .env
set +o allexport

# ....Helper function..............................................................................
# import shell functions from utilities library
N2ST_PATH=${N2ST_PATH:-"${NABO_PATH}/build_system/utilities/norlab-shell-script-tools"}
source "${N2ST_PATH}/import_norlab_shell_script_tools_lib.bash"

# ====Begin========================================================================================
print_formated_script_header 'nabo_execute_nabo_unittest.bash' ':'

cd "${NBS_LIB_INSTALL_PATH}/${NBS_REPOSITORY_NAME}/build"

if [[ ${IS_TEAMCITY_RUN} == true ]] || [[ ${TEAMCITY_VERSION} ]]; then
  echo -e "##teamcity[testSuiteStarted name='gtest']"
  echo -e "##teamcity[testStarted name='gtest' captureStandardOutput='true']"
else
  print_msg "Starting Libnabo GoogleTest unit-test"
fi

# .................................................................................................
sudo chmod +x test
make test
NABO_TEST_EXIT_CODE=$?
# .................................................................................................

SUCCESS_MSG="Libnabo GoogleTest unit-test completed successfully"
FAILURE_MSG="Libnabo GoogleTest unit-test completed with error"

if [[ ${IS_TEAMCITY_RUN} == true ]] || [[ ${TEAMCITY_VERSION} ]]; then
  echo -e "##teamcity[testFinished name='gtest']"

  # Report message to build log
  if [[ ${NABO_TEST_EXIT_CODE} == 0 ]]; then
    echo -e "##teamcity[message text='${MSG_BASE_TEAMCITY} ${SUCCESS_MSG}' status='NORMAL']"
  else
    echo -e "##teamcity[message text='${MSG_BASE_TEAMCITY} ${FAILURE_MSG}' errorDetails='$NABO_TEST_EXIT_CODE' status='ERROR']" 1>&2
  fi

  echo -e "##teamcity[testSuiteFinished name='gtest']"
else

  if [[ ${NABO_TEST_EXIT_CODE} == 0 ]]; then
    print_msg_done "${SUCCESS_MSG}"
  else
    print_msg_error "${FAILURE_MSG}"
  fi

fi

print_formated_script_footer 'nabo_execute_nabo_unittest.bash' ':'
# ====Teardown=====================================================================================
cd "${TMP_CWD}"
exit $NABO_TEST_EXIT_CODE
