libnabo is a fast K Nearest Neighbour library for low-dimensional spaces.
It provides a clean, legacy-free, scalar-type–agnostic API thanks to C++ templates.
Its current CPU implementation is strongly inspired by [ANN], but with more compact data types.
On the average, libnabo is 20% faster than [ANN].

libnabo depends on [Eigen], a modern C++ matrix and linear-algebra library.
libnabo is being developed by Stéphane Magnenat as part of his work at [ASL-ETH](http://www.asl.ethz.ch).


Compilation
-----------

libnabo uses [CMake] as build system.
Just create a directory, go inside it and type:

	cmake LIBNABO_SRC_DIR
    
where `LIBNABO_SRC_DIR` is the top-level directory of libnabo's sources.
If [Eigen] is not installed system wide, you might have to tell [CMake] where to find it.
Please read the [CMake documentation].


Usage
-----

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


Unit testing
------------

The distribution of libnabo integrates a unit test module, based on CTest.
Just type:

	make test
   
...in the build directory to run the tests.
Their outputs are available in the `Testing` directory.
These consist of validation and benchmarking tests.
If [ANN] or [FLANN] are detected when compiling libnabo, `make test` will also perform comparative benchmarks.


Bug reporting
-------------

Please use [github's issue tracker](http://github.com/ethz-asl/libnabo/issues) to report bugs.

License
-------

libnabo is released under a permissive BSD license.


[ANN]: http://www.cs.umd.edu/~mount/ANN
[FLANN]: http://www.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN
[CMake]: http://www.cmake.org
[CMake documentation]: http://www.cmake.org/cmake/help/cmake2.6docs.html
[Eigen]: http://eigen.tuxfamily.org