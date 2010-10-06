libnabo is a fast K Nearest Neighbour library for low-dimensional spaces.
It provides a clean, legacy-free, scalar-type–agnostic API thanks to C++ templates.
Its current CPU implementation is strongly inspired by [ANN], but with more compact data types.
On the average, libnabo is 20% faster than ANN.

libnabo depends on [Eigen], a modern C++ matrix and linear-algebra library.

Usage
-----

libnabo is easy to use. For example, assuming that you are working with floats and that you have a point set `M` and a query point `q`, you can find the indexes `n` of the K nearest neighbours of `q` in `M`:

	#include "nabo/nabo.h"
	using namespace Nabo;
	using namespace Eigen;
	...
	NNSearchF* nns = NNSearchF::createKDTreeLinearHeap(M);
	const int K(5);
	VectorXi n = nns->knn(q, K);

In this example, `M` is an [Eigen] matrix (column major, float) and `q` is an [Eigen] vector (float).
The result `n` is an [Eigen] vector of indices of columns of `M`.
See `example/trivial.cpp` for a compilable version of this example.

More information available soon.

Compilation
-----------

libnabo uses [CMake] as build system.
Just create a directory, go inside it and type:

	cmake LIBNABO_SRC_DIR
    
where `LIBNABO_SRC_DIR` is the top-level directory of libnabo's sources.

Testing
-------

The distribution of libnabo integrates a unit test module, based on CTest.
Just type:

	make test
   
...in the build directory to run the tests.
Their outputs are available in the `Testing` directory.
These consist of validation and benchmarking tests.
If [ANN] is detected when compiling libnabo, `make test` will also perform comparative benchmarking.

[ANN]: http://www.cs.umd.edu/~mount/ANN
[CMake]: http://www.cmake.org
[Eigen]: http://eigen.tuxfamily.org

