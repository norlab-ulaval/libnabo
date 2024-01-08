<div align="center">

[//]: # (====GitHub badges========================================================================)

<img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/norlab-ulaval/libnabo">
<img alt="GitHub forks" src="https://img.shields.io/github/forks/norlab-ulaval/libnabo">
<img alt="GitHub License" src="https://img.shields.io/github/license/norlab-ulaval/libnabo">
<img alt="GitHub release (with filter)" src="https://img.shields.io/github/v/release/norlab-ulaval/libnabo">
<img src="https://img.shields.io/static/v1?label=JetBrains TeamCity&message=CI/CD&color=green?style=plastic&logo=teamcity" />
<br>
<br>


[//]: # ( ==== Description =========================================== )
**libnabo is a fast K Nearest Neighbour library for low-dimensional spaces.<br>
It provides a clean, legacy-free, scalar-type–agnostic API thanks to C++ templates.**
<br>
Its current CPU implementation is strongly inspired by [ANN], but with more compact data types.<br>
On the average, libnabo is 5% to 20% faster than [ANN].
<br>
<br>

</div>

libnabo depends on [Eigen], a modern C++ matrix and linear-algebra library.
libnabo works with either version 2 or 3 of Eigen.
libnabo also optionally depends on [Boost], a C++ general library, for Python bindings.

libnabo was developed by [Stéphane Magnenat](http://stephane.magnenat.net) as part of his work at [ASL-ETH](http://www.asl.ethz.ch) and is now maintained by [Simon Lynen](https://github.com/simonlynen).

If you are interested in a pure-[Rust](https://www.rust-lang.org/) version, check [that repository](https://github.com/enlightware/nabo-rs) out.

---

[//]: # (====Supported OS and aarch===============================================================)

### Supported OS And Architecture
libnabo is tested on our build system under the following architecture and OS:

- x86 and arm64/v8
- Ubuntu bionic (18.04) focal (20.04) and jammy (22.04)

Note:

- libnabo reportedly works on MacOs OsX (latest) and Windows (latest)

---

[//]: # (====Release note=========================================================================)

### ★ Version `1.1.0` Release Note 

This release of _libnabo_ introduces the integration
of [norlab-build-system (NBS)](https://github.com/norlab-ulaval/norlab-build-system) as a _git
submodule_ for codebase development and testing.

Execute the following to clone the repository with its submodule:

```shell
git clone --recurse-submodules https://github.com/norlab-ulaval/libnabo.git
```

If _libnabo_ was previously cloned, execute the following to fetch its new submodule

```shell
git submodule update --remote --recursive --init
```

### ★ Contributing Instructions

See [contributing_instructions.md](contributing/contributing_instructions.md)
for instructions related to bug reporting, code contribution and for setting up
the `libnabo-build-system`
on your workstation to speed up your local development workflow.


---


Download
========

Ubuntu builds are available on my PPA at: https://launchpad.net/~stephane.magnenat
They provide a package with the shared library, another with the development headers and a third with the documentation.

The source code is available from github, you can clone the git tree by doing:

```shell
git clone --recurse-submodules https://github.com/norlab-ulaval/libnabo.git
```


Docker images 
========
Run the following commands to pull and run libnabo in a docker container 
```shell
docker pull norlabulaval/libnabo:latest-ubuntu-focal

docker run -it --rm norlabulaval/libnabo:latest-ubuntu-focal
```
See available [libnabo image tags](https://hub.docker.com/repository/docker/norlabulaval/libnabo/) on dockerhub.

To install docker related dependencies on ubuntu, execute the following

```shell
cd ./build_system/nabo_utility_script

# Execute docker tools install script i.e. docker daemon, docker compose, docker buildx
bash nabo_install_docker_tools.bash
```

Compilation
===========

libnabo uses [CMake] as build system.
The complete compilation process depends on the system you are using (Linux, Mac OS X or Windows).
You will find a nice introductory tutorial in [this video](http://www.youtube.com/watch?v=CLvZTyji_Uw).

For conveniences, you can use the provided installer script for ubuntu
```shell
bash libnabo_dependencies_installer.bash

# Use the --help flag to see the list of optional flag
bash libnabo_installer.bash [<optional flag>]
```

Prerequisites
-------------

If your operating system does not provide it, you must get [Eigen], and [Boost] if you want to build the Python bindings.
[Eigen] only needs to be downloaded and extracted.

Compilation options
-------------------

libnabo provides the following compilation options, available through [CMake]:

 * `SHARED_LIBS` (boolean, default: `false`): if `true`, build a shared library, otherwise build a static library

You can specify them with a command-line tool, `ccmake`, or with a graphical tool, `cmake-gui`.
Please read the [CMake documentation] for more information.

Quick compilation and installation under Unix
---------------------------------------------

Under Unix, assuming that [Eigen] and [Boost] are installed system-wide, you can compile (with optimisation and debug information) and install libnabo in `/usr/local` with the following commands run in the top-level directory of libnabo's sources:

	SRC_DIR=`pwd`
	BUILD_DIR=${SRC_DIR}/build
	mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ${SRC_DIR}
	# if Eigen or Boost are not available system-wide, run at that point:
	#   cmake-gui .
	# cmake-gui allows you to tell the location of Eigen or Boost
	make
	sudo make install

These lines will compile libnabo in a `build` sub-directory and therefore keep your source tree clean.
Note that you could compile libnabo anywhere you have write access, such as in `/tmp/libnabo`.
This out-of-source build is a nice feature of [CMake].
If [Eigen] or [Boost] are not installed system-wide, you might have to tell [CMake] where to find them (using `ccmake` or `cmake-gui`).

You can generate the documentation by typing:

	make doc

Usage
=====

libnabo is easy to use. For example, assuming that you are working with floats and that you have a point set `M` and a query point `q`, you can find the `K` nearest neighbours of `q` in `M`:

	#include "nabo/nabo.h"
	using namespace Nabo;
	using namespace Eigen;
	...
	NNSearchF* nns = NNSearchF::createKDTreeLinearHeap(M);

	const int K = 5;
	VectorXi indices(K);
	VectorXf dists2(K);

	nns->knn(q, indices, dists2, K);

In this example, `M` is an [Eigen] (refering to the software, not to the math) matrix (column major, float) and `q` is an [Eigen] vector (float). Note that `M` **must stay alive** throughout the use of libnabo, otherwise the results of `knn` are undefined.
The results `indices` and `dists2` are [Eigen] vectors of indices and squared distances refering to the columns of `M`.
See `examples/trivial.cpp` for a compilable version of this example, and `examples/usage.cpp` for a slightly more complex example involving multi-point queries.

Running `make doc` in your build directory will generate a browsable documentation in `doc/html`.
The main page `doc/html/index.html` contains a detailed overview of the usage of libnabo.

You can find a complete CMake integration
example in [examples/libnabo-cmake-example](examples/libnabo-cmake-example) to
see how to look for, and link against this library.

Python bindings
===============

libnabo includes python bindings that are compiled if Python is available.
The resulting module is called pynabo, you can see an example in `python/test.py`.
You can find more information in the docstring-based documentation:

	python -c "import pynabo; help(pynabo.NearestNeighbourSearch)"

Building
--------

The Python bindings can be generated for Python 2 or Python 3.
To specify the version of the interpreter to use when building the bindings, set the `PYTHON_VERSION_MAJOR` and `PYTHON_VERSION_MINOR` variables.
For example if you have both Python 2.7 and 3.5 installed, you could ask CMake to generate Python 3 bindings by using the following command.

    cmake -DPYTHON_VERSION_MAJOR=3 -DPYTHON_VERSION_MINOR=5 ..

On Debian-based distributions you may also need the `-DPYTHON_DEB_INSTALL_TARGET` option enabled.

Unit testing
============

The distribution of libnabo integrates a unit test module, based on CTest.
Just type:

	make test

...in the build directory to run the tests.
Their outputs are available in the `Testing` directory.
These consist of validation and benchmarking tests.
If [ANN] or [FLANN] are detected when compiling libnabo, `make test` will also perform comparative benchmarks.

Citing libnabo
==============

If you use libnabo in the academic context, please cite this paper that evaluates its performances in the contex of ICP:

	@article{elsebergcomparison,
		title={Comparison of nearest-neighbor-search strategies and implementations for efficient shape registration},
		author={Elseberg, J. and Magnenat, S. and Siegwart, R. and N{\"u}chter, A.},
		journal={Journal of Software Engineering for Robotics (JOSER)},
		pages={2--12},
		volume={3},
		number={1},
		year={2012},
		issn={2035-3928}
	}



License
=======

libnabo is released under a permissive BSD license.


[ANN]: http://www.cs.umd.edu/~mount/ANN
[FLANN]: http://www.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN
[CMake]: http://www.cmake.org
[CMake documentation]: http://www.cmake.org/cmake/help/cmake2.6docs.html
[Eigen]: http://eigen.tuxfamily.org
[Boost]: http://www.boost.org
