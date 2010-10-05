libnabo is a fast K Nearset Neighbor library for low-dimensional spaces.
It provides a clean, legacy-free, scalar-type–agnostic API thanks to C++ templates.
Its current CPU implementation is strongly inspired by ANN [1], but with more compact data types.
On the average, libnabo is 20% faster than ANN.

Usage
-----

TODO

Compilation
-----------

libnabo uses CMake [2] as build system.
Just create a directory, go inside it and do:

	cmake LIBNABO_SRC_DIR
    
where `LIBNABO_SRC_DIR` is the top-level directory of libnabo's sources.

Testing
-------

The distribution of libnabo integrates a unit test module, based on CTest.
Just type:

	make test
   
...in the build directory to run the tests.
Their outputs are available in the `Testing` directory.

Benchmarking
------------

If ANN [1] is detected when compiling libnabo, benchmarking will be available

[1]: http://www.cs.umd.edu/~mount/ANN/
[2]: http://www.cmake.org/

