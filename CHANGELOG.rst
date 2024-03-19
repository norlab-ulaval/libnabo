^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package libnabo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Forthcoming
-----------

* Update README.md
* Set minimum required CMake version to 3.10.2
* (origin/dev-feat-integrate-norlab-build-system-library-and-implement-linabo-build-system-NMO-405) style: renamed `.*test_compilation.*` services and target to `.*integration_test.*` for clarity.
* doc(readme): mute ubuntu jammy from supported OS version until libpointmatcher build pass on jammy
* ci: drop ubuntu jammy from build matrix
* doc: add missing norlab light logo
* doc: add missing norlab light logo
* doc: update header with dynamic logo, fix hyperlink, add dockerhub badge and fix relative link
* doc(pull_request_template.md): fix typo
* doc(readme): add norlab logo to header and improve intro
* test: add release crawler script bats test
* build: add release crawler script and .env file
* doc(readme): update header
* build: update repo version format
* doc: clean pr comment
* feat: integrate norlab-build-system library and implement linabo-build-system [NMO-405]
* build: Added norlab-shell-script-tools submodule to repository [NMO-405]
* build: Added norlab-build-system submodule to repository [NMO-405]
* chore: added a code owner designation file, a pull request template and a conventional commit reference file
* Update README.md
* Update README.md with cmake guide
* Use CMakePackageConfigHelpers to generate config files. Fill libnabo_INCLUDE_DIRS and libnabo_LIBRARIES variables.
* Fix Cmake Syntax invoked by Nabo's version regex. Updated min cmake version to 3.8
* more robust extraction of version
* Remove dependency on grep + Improve findstr command
* Use findstr instead of grep on Windows
* Date range and catch-all for authors
* Upgrade all syntax to package format 3
* format "3" specifier for condition flag
* Keep catkin for ROS1
* Create LICENSE file based on BSD license as per package.xml
* catkin not required for pure cmake packages
* Added link to Rust version.
* Replace TABS with SPACES
* Automaticaly find eigen3
* Fix install-space include directories
* fix missing eigen headers in python module
* add complete cmake exported-targets usage example
* option to disable doxygen; cmake clean up
* better handling of STATIC build option
* modernize cmake to have exported targets
* Make library git submodule-friendly

1.0.7 (2019-02-07)
------------------
* Disabled cmake compile tests by default and on compilers that do not support them (#95)
* Fix Python 2 bindings support in CMake scripts (#90)
* Port libnabo to c++11 (#89)
* Remove register keyword in index_heap.h (#88)
* Fix compilation warning for MSVC (#85)
* Assert template type for invalid setters (#80)
* Return numerically maximal index (unsigned) or -1 (signed) for no match case (#79)
* Add generate step for ${PROJECT_BINARY_DIR}/libnaboConfig.cmake (#76)
* Removed compiler-specific flags for compilers that do not support them (#74)
* Added cmake_policy(SET CMP0054 NEW) (#73)
* Output compiler message with compile test fatal error in cmake (#72)
* Removed erroneous commas from test/CMakeLists.txt (#71)
* Removed fatal "," suffix from FATAL_ERROR in CMakeLists.txt (#70)
* Fixed regression concerning installed libnaboConfig.cmake (#65)
* Fix/relax compiler requirements (#63)
* Removed hard dependency on the doc target (#62)
* Install any.hpp (#61)
* Remove boost::any and boost:format dependencies (#59)
* Port the python bindings to python3 (#57)
* Added cmake switch to disable usage of OpenMP (#53)
* Zero copy for Eigen::Matrix3XT and Eigen::Map<const Eigen::Matrix3XT> (#43)
* Fix warnings and switch on Wextra (#42)
* Disallow instantiation with non dynamic matrices (#41)
* Update README.md
* Removed all code dealing with libnaboTargets.cmake (#32)
* Got rid of unused locally defined typedefs (#27)
* Contributors: David Landry, Hannes Sommer, Simon Lynen, Simon-Pierre Deschênes, Stéphane Magnenat, cezheng, ffurrer, magehrig, renning22, sandsmark, taketwo, tcies

1.0.6 (2015-03-05)
------------------
* Reset point indices of results with distances exceeding threshold (#23, #24)
* Fine tune the find_package() capability and add uninstall target (#22)
* Fixed compiler warning (#18)
* Added OpenMP support (#20, #21)
* Build type tuning (#19)
* Fix: terminal comma in enum requires C++11
* Fix UBSAN error calculating maxNodeCount (#16, #17)
* Fixed tiny (yet significant) error in the Python doc strings (#15)
* Compile static lib with PIC (#14)
* Contributors: Francois Pomerleau, François Pomerleau, Gregory Hitz, Gregory Jefferis, Simon Lynen, Stéphane Magnenat

1.0.5 (2014-06-12)
------------------
* Added configure scripts for full catkinization
* Catkinization of libnabo (following REP136)
* Update README.md
  Added Simon as the maintainer.
* [test] use CLOCK_PROF for NetBSD build
* Fixed CppCheck warning.
  Fix broken install when doxygen is not found
* Fix cmake stylistic issue
* Make python install respect custom CMAKE_INSTALL_PREFIX
* Fix broken install when doxygen is not found
* Contributors: Chris Foster, Francis Colas, Paul Furgale, Pierrick Koch, Stéphane Magnenat, fcolas

1.0.4 (2013-09-03)
------------------
* Updated Debian version number.
* Added check for invalid optionFlags values.
* Fixed compilation on OS X, reduced code duplication in tests by moving timers into helpers.h
* Contributors: Stéphane Magnenat

1.0.3 (2013-08-19)
------------------
* Prevent requesting more points than available in the cloud, prevent the use of empty clouds, bumped version number.
* Added test for grep.
* Worked around issue `#3 <https://github.com/ethz-asl/libnabo/issues/3>`_ on Windows.
* Updated documentation.
* Contributors: Stéphane Magnenat

1.0.2 (2013-02-20)
------------------
* Updated version number.
* Added const
* Removed useless optimisation, cleaned-up.
* faster tree heap when less neighbors than requested
* Hopefully fixed detection of missmatched python versions.
* Fixed cloud size check for clouds with billions of points.
* Contributors: Francis Colas, Stéphane Magnenat

1.0.1 (2012-11-06)
------------------
* Bumped version in doc.
* Updated documentation to reflect the disabling of OpenCL.
* Fixed bug in version define.
* Typo.
* Added note about nabo citation.
* Fixed the semantic of epsilon to match the documentation's and ANN's ones.
* Added ability to pass one radius per point.
* Added custom command for python target.
* Reverted buggy change
* Fixed old python_add_module
* fixed python 2.6 detection on debian.
* Fixed debian installation.
* Added debian-specific install target.
* Fixed doc and python link command
* Added docstring to python bindings.
* Added test and improved CMakeLists.txt
* Added debian install for python
* Renammed to pynabo, added python install script.
* Experimental python bindings working.
* When number of point is smaller than bucket size, create a single-bucket tree.
* Cleaned-up makefile.
* Fixed description.
* Merge branch 'master' of github.com:ethz-asl/libnabo
* Added dbg package.
* Contributors: Stéphane Magnenat

1.0.0 (2011-10-19)
------------------
* Fixed naming convention.
* Removed dbg msg.
* Fixed doc.
* Updated doc, fixed debian build.
* Updated README.
* Fixed debian compilation.
* Added note about download.
* Fixed bug in control file.
* Separated doc package from dev package.
* Cleaned-up debian build chain.
* Added debian package.
* Added link to online doc.
* Added major version in library name.
* Fixed documentation.
* Improved Makefile and documentation.
* Added bench to select bucket size.
* Updated README.
* Fixed doxygen warning.
* Minor changes.
* Been kind to Francis and in example compile in a build subdirectory.
* Minor fix
* Improved documentation.
* Search for Eigen in ROS diamondback by default.
* Updated (c) date.
* Added using directive for boost.
* Added win32 compatibility (thanks Alessio Placitelli)
* Fixed bug when dimension was not passed.
* Added const to knn search, bumped version number.
* Removed duplicated comment.
* Improved documentation.
* Added additional search parameters to specify bucketSize for CPU kd-trees.
* Optimized memory structure for CPU-basde kd-tree.
* Added buckets.
* Added radius search.
* Fixed test case when CL is disabled. Improved verbose output of configuration.
* Cleaned-up OpenGL API, marked it as unstable.
* Fixed compilation of OpenCL part. Added high-res timer for benches when available.
* Search for eigen (3) not explicitely eigen 2.
* Added Eigen3 compatibility.
* Result-file header now has the right number of columns.
* Added statistics infrastructure.
* Added caching to OpenCL
* Removed arbitrary constant before method.
* Added missing files.
* Added epsilon test.
* Added link to FLANN
* Fixed typo
* Fixed link
* Added virtual destructor to NNS interface to prevent memory leak in children.
* Fixed clang compilation.
* Fixed extraction of version
* Added new method for GPU-based kd tree.
* removed dependency on C++0x
* Updated to latest draft of C++0x
* Fixed compilation when OpenCL is not present
* Updated doc
* Merge branch 'master' of github.com:stephanemagnenat/libnabo
* Fixed implementation to fit new API.
* Changed API. Implementation broken.
* Fix compilation with undefined HAVE_OPENCL
* changes names of variables to avoid overlaying.
* Added multiple query per run.
* Fixed uninitialised memory.
* Fixed buffer handling for OpenCL, there seems to be still a bug with memory.
* Improved OpenCL infrastructure.
* Added back GPU
* OpenCL KDtree now working.
* Fixed adresse of node array.
* OpenCL kernel for NNS compiles.
* OpenCL glue now works to the point of reporting compilation errors in the source code.
* Written OpenCL kernel for knn search, glue is still needed.
* Added infrastructure for OpenCL support.
* Added flann comparison
* improved diff to ANN
* improved doc.
* Added documentation to source code
* Improved README.
* Added more complex example.
* Added license
* Improved README.
* Added example
* Restructured library.
* Improved readme.
* added initial readme.
* Use index instead of values for temporary vector to create nodes, results in a faster creation.
* libnabo now always faster than ANN.
* Cleaned-up bench infrastructure, now it is possible to do more than one time each bench.
* Added reentrant statistics, depends on C++0x.
* Added ref to points in dist function, equals perfs of ANN.
* Fixed KDTree.
* Added explicit bound version of KDTree, ANN style.
* Added unbalanced tree.
* ANN bench now has both search and pri-search
* Improved bench API.
* Added option for cell balancing.
* Prevent overflow in stats.
* Improved benchmarking.
* Added pt in leave option.
* Added stack-based KNN on our structure, same perf as priority_queue... still 2x worst than ANN, memory-bounded?
* Added API to match several points at once.
* Added bench, comparison with ANN
* Added large test.
* Improved tests
* Added unit tests.
* Fixe includes for Lucid's version of Eigen lib.
* use better dist
* fixed bug
* Restructured project.
* added missing files
* Refactored API.
* Improved performance of search in kdtree.
* Renamed lib, should help compilation with old cmakes.
* Fixed arbitrary dimensions.
* Added search in kdtree
* Contributors: Francois Pomerleau, Martin Voelkle, Stéphane Magnenat
