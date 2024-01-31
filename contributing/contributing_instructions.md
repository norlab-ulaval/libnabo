# Contributing to _libnabo_

## Bug Reporting

Please use our [github's issue tracker](http://github.com/norlab-ulaval/libnabo/issues) to
report bugs.

## Code Contributions

Libnabo codebase now
integrate [norlab-build-system (NBS)](https://github.com/norlab-ulaval/norlab-build-system)
and [norlab-shell-script-tools (N2ST)](https://github.com/norlab-ulaval/norlab-shell-script-tools).
`NBS` is a build-infrastructure-agnostic build system custom-made to meet our needs in robotic
software engineering at NorLab and `N2ST` is a library of shell script functions as well as a shell
testing tools leveraging _**bats-core**_ and _**docker**_ .
`N2ST` purpose is to speed up shell script development and improve reliability.

`NBS` is deployed on our [TeamCity](https://www.jetbrains.com/teamcity/) continuous
integration/deployment server and oversees protected branches of
the [libnabo](https://github.com/norlab-ulaval/libnabo) GitHub repository:

- The `develop` branch can only be merged through a pull-request from any `<feature>` branches. Any
  contributor can submit a pull request to the `develop` branch;
- the `release` branch is a revision and preparation branch where we can freeze the codebase in a
  given state without stalling to `develop` branch progression;
- The `master` branch can only be merged through a pull-request from the `release` branch. Only
  repository admin can submit a PR to the `master` branch.

In any cases, submitting a pull request to `develop` or `master` will trigger a build/test
configuration on our build system and the pull request will be granted if the build/test run
succeed.

**Current build matrix:**
`[latest] x [x86, arm64] x [ubuntu] x [bionic, focal] x [Release, RelWithDebInfo, MinSizeRel]`

### Development Workflow

To speed up the development process, you can run the build system localy on your workstation and
have access to stacktrace and build log.
It support multi-OS and multi-architecture through docker container.

#### Install _libnabo-build-system_ Dependencies

```shell
cd <path/to/libnabo>

# If libnabo is already cloned, fetch the NBS and N2ST submodule 
git submodule update --remote --recursive --init

cd ./build_system/nabo_utility_script

# Execute docker tools install script i.e. docker daemon, docker compose, docker buildx 
bash nabo_install_docker_tools.bash

# Configure a multi-architecture docker builder
bash nabo_create_multiarch_docker_builder.bash
```

#### _libnabo_ Development › To Execute Build/Test Step Locally

```shell
cd <path/to/libnabo>/build_system

# Run the build matrix as specified in ".env.build_matrix.libnabo" 
#   on native architecture using "ci_PR" service 
bash nabo_crawl_libnabo_build_matrix.bash --fail-fast -- build ci_PR

# Run a specific case using build flags with multi-architecture 
# virtualization using "ci_PR_amd64" and "ci_PR_arm64v8" services 
bash nabo_crawl_libnabo_build_matrix.bash \
            --repository-version-build-matrix-override latest \
            --os-name-build-matrix-override ubuntu \
            --cmake-build-type-build-matrix-override RelWithDebInfo \
            --ubuntu-version-build-matrix-override focal \
            --fail-fast \
            -- build ci_PR_amd64 ci_PR_arm64v8

# Read the help for details
bash nabo_crawl_libnabo_build_matrix.bash --help
```

Note: To assess the state of the codebase, even for cases that are known the break the build,
execute `nabo_crawl_libnabo_build_matrix.bleeding.bash` with build
matrix `.env.build_matrix.libnabo.bleeding`.
The stable build matrix used for release is `.env.build_matrix.libnabo`.

#### Build System Development

```shell
cd <path/to/libnabo>/build_system/tests/
 
# To execute docker dryrun and configuration tests
bash run_all_docker_dryrun_and_config_tests.bash

# To execute shell script tests
bash run_bats_core_test_in_n2st.bash

# To spin a container in interactive mode with the codebase cloned but not compiled  
cd ./tests_docker_interactive/
bash build_and_run_IamBuildSystemTester.bash bash
```

#### Build System Notes

- `nabo_crawl_dependencies_build_matrix.bash` execute the build matrix for the libnabo
  dependencies.
  It's not required to build them locally as they are pre-build by our TeamCity server periodically
  push to dockerhub.
  When executing `nabo_crawl_libnabo_build_matrix.bash`, the `libnabo-dependencies`
  docker images are pull and used as base image for the `libnabo-[ci_PR_test|release]`
  images.
- About `libnabo/.github/workflow/` vs `libnabo/build_system/` logic: Those are
  separate build logic.
  `.github/workflow/` was community contributed and as the responsibilities of building
  python-binding and pushing packages.
  For this reason, it run a one-dimension build matrix: multiple python version, single OS version,
  single arch (x86) and
  single compile flag which GitHub action computing resources can handle just fine.

## Commit Messages

This is optional for now but will eventually move our release workflow to semantic-versioning.
See [Commit Message References](commit_msg_reference.md) for details.

## Note For Repository Admins

### About Release Branch And Pull Request To Master Branch

- Only repository admins have the privilege to `push/merge` on the default branch (ie: `master`)
  and the `release` branch.
- Keep PR in `draft` mode until all the release reviewers are ready to push the release.
- Once a PR from `release` -> `master` branch is created (not in draft mode),
    - it triggers the _build-system_ test
    - (in-progress) and it triggers the _semantic release automation_
