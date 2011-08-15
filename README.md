libnabo is a fast K Nearest Neighbour library for low-dimensional spaces.
It provides a clean, legacy-free, scalar-type–agnostic API thanks to C++ templates.
Its current CPU implementation is strongly inspired by [ANN], but with more compact data types.
On the average, libnabo is 5% to 20% faster than [ANN].

libnabo depends on [Eigen], a modern C++ matrix and linear-algebra library.
libnabo works with either version 2 or 3 of Eigen.
libnabo is being developed by [Stéphane Magnenat](http://stephane.magnenat.net) as part of his work at [ASL-ETH](http://www.asl.ethz.ch).


Compilation
===========

libnabo uses [CMake] as build system.
The complete compilation process depends on the system you are using (Linux, Mac OS X or Windows).
You will find a nice introductory tutorial in [this video](http://www.youtube.com/watch?v=CLvZTyji_Uw).

Quick compilation and installation under Unix
---------------------------------------------

Under Unix, assuming that [Eigen] is installed system-wide, you can compile (with optimisation and debug information) and install libnabo in `/usr/local` with the following commands run in the top-level directory of libnabo's sources:

	SRC_DIR=`pwd`
	BUILD_DIR=${SRC_DIR}/build
	mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ${SRC_DIR}
	make
	sudo make install

These lines will compile libnabo in a `build` sub-directory and therefore keep your source tree clean.
Note that you could compile libnabo anywhere you have write access, such as in `/tmp/libnabo`.
This out-of-source build is a nice feature of [CMake] under Unixes.

If [Eigen] is not installed system-wide, you might have to tell [CMake] where to find it.
You can do this with a command-line tool, `ccmake`, or with a graphical tool, `cmake-gui`.
Please read the [CMake documentation] for more information.


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

In this example, `M` is an [Eigen] (refering to the software, not to the math) matrix (column major, float) and `q` is an [Eigen] vector (float).
The results `indices` and `dists2` are [Eigen] vectors of indices and squared distances refering to the columns of `M`.
See `examples/trivial.cpp` for a compilable version of this example, and `examples/usage.cpp` for a slightly more complex example involving multi-point queries.

Running `make doc` in your build directory will generate a browsable documentation in `doc/html`.
The main page `doc/html/index.html` contains a detailed overview of the usage of libnabo.


Unit testing
============

The distribution of libnabo integrates a unit test module, based on CTest.
Just type:

	make test
   
...in the build directory to run the tests.
Their outputs are available in the `Testing` directory.
These consist of validation and benchmarking tests.
If [ANN] or [FLANN] are detected when compiling libnabo, `make test` will also perform comparative benchmarks.


Bug reporting
=============

Please use [github's issue tracker](http://github.com/ethz-asl/libnabo/issues) to report bugs.


License
=======

libnabo is released under a permissive BSD license.


[ANN]: http://www.cs.umd.edu/~mount/ANN
[FLANN]: http://www.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN
[CMake]: http://www.cmake.org
[CMake documentation]: http://www.cmake.org/cmake/help/cmake2.6docs.html
[Eigen]: http://eigen.tuxfamily.org