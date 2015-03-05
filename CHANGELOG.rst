^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package libnabo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
