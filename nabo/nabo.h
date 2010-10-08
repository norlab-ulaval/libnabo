/*

Copyright (c) 2010, Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the author at <stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef __NABO_H
#define __NABO_H

#include "Eigen/Core"
#include "Eigen/Array"
#include <vector>
#include <cstdatomic>

/*! 
	\file nabo.h
	\brief public interface
	\ingroup public
*/

/*!
\mainpage libnabo

from http://github.com/stephanemagnenat/libnabo by Stéphane Magnenat (http://stephane.magnenat.net),
ASL-ETHZ, Switzerland (http://www.asl.ethz.ch)

libnabo is a fast K Nearest Neighbour library for low-dimensional spaces.
It provides a clean, legacy-free, scalar-type–agnostic API thanks to C++ templates.
Its current CPU implementation is strongly inspired by \ref ANN, but with more compact data types.
On the average, libnabo is 20% faster than \ref ANN.

libnabo depends on \ref Eigen, a modern C++ matrix and linear-algebra library.

\section Compilation

libnabo uses \ref CMake as build system.
Just create a directory, go inside it and type:

\code
cmake LIBNABO_SRC_DIR
\endcode

where \c LIBNABO_SRC_DIR is the top-level directory of libnabo's sources.
If \ref Eigen is not installed system wide, you might have to tell \ref CMake where to find it.
Please read the <a href="http://www.cmake.org/cmake/help/cmake2.6docs.html">CMake documentation</a>.

\section Usage

libnabo is easy to use. For example, assuming that you are working with floats and that you have a point set \c M and a query point \c q, you can find the indices \c n of the K nearest neighbours of \c q in \c M :

\include trivial.cpp

In this example, \c M is an \ref Eigen (refering to the software, not to the math) matrix (column major, float) and \c q is an \ref Eigen vector (float).
The result \c n is an \ref Eigen vector of indices refering to the columns of \c M.

Here is a slightly more complex example:

\include usage.cpp


\section UnitTesting Unit testing

The distribution of libnabo integrates a unit test module, based on CTest.
Just type:

\code
make test
\endcode
   
...in the build directory to run the tests.
Their outputs are available in the \c Testing directory.
These consist of validation and benchmarking tests.
If \ref ANN is detected when compiling libnabo, \c make \c test will also perform comparative benchmarks.


\section BugReporting Bug reporting

Please use <a href="http://github.com/stephanemagnenat/libnabo/issues">github's issue tracker</a> to report bugs.

\section License

libnabo is released under a permissive BSD license.

\section Faq

\subsection ANN

libnabo differs from \ref ANN on the following points:

* API
- templates for scalar types
- self-match option as execution-time (instead of compile-time) parameter
- Eigen library [2] for vector and matrixes

* limitations
- currently no radius search
- currently only euclidean distance
- currently only one-point buckets

* implementation
- optional O(log(n)) tree heap instead of O(n) vector heap
- compact memory representation, one memory allocation for all nodes
- implicit reference to left child (always next node in array)
- do not store bounds in nodes (that is, I do it like in your article)

* performances
- about 20% faster than ANN (both -O3 -NDEBUG)
- clearly memory-bound, neither OpenMP nor boost::thread improve performances 

\section References

\li \anchor Eigen Eigen: http://eigen.tuxfamily.org
\li \anchor ANN ANN: http://www.cs.umd.edu/~mount/ANN
\li \anchor CMake CMake: http://www.cmake.org

*/

//! Namespace for Nabo
namespace Nabo
{
	//! \defgroup public public interface 
	//@{
	
	//! version of the Nabo library as string
	#define NABO_VERSION "0.9.0"
	//! version of the Nabo library as an int
	#define NABO_VERSION_INT "9000"
	
	//! Nearest neighbour search interface, templatized on scalar type
	template<typename T>
	struct NearestNeighbourSearch
	{
		//! an Eigen vector of type T, to hold the coordinates of a point
		typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vector; 
		//! a column-major Eigen matrix in which each column is a point; this matrix has dim rows
		typedef typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;
		//! an index to a Vector or a Matrix, for refering to data points
		typedef int Index;
		//! a vector of indices to data points
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		//! a matrix of indices to data points
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> IndexMatrix;
		
		//! the reference to the data-point cloud, which must remain valid during the lifetime of the NearestNeighbourSearch object
		const Matrix& cloud;
		//! the dimensionality of the data-point cloud
		const Index dim;
		//! the low bound of the search space (axis-aligned bounding box)
		const Vector minBound;
		//! the high bound of the search space (axis-aligned bounding box)
		const Vector maxBound;
		
		//! statistics for search
		struct Statistics
		{
			//! create zero-initialised statistics
			Statistics():lastQueryVisitCount(0),totalVisitCount(0) {}
			std::atomic_uint lastQueryVisitCount; //!< number of visits during the last query
			std::atomic_uint totalVisitCount; //!< total number of visits
		};
		
		//! type of search
		enum SearchType
		{
			BRUTE_FORCE = 0, //!< brute force, check distance to every point in the data
			KDTREE_LINEAR_HEAP = 1, //!< kd-tree with linear heap, good for small k (~up to 30)
			KDTREE_TREE_HEAP = 2, //!< kd-tree with tree heap, good for large k (~from 30)
			SEARCH_TYPE_COUNT //!< number of search types
		};
		
		//! option
		enum SearchOptionFlags
		{
			ALLOW_SELF_MATCH = 1, //!< allows the return of the same point as the query, if this point is in the data cloud; forbidden by default
			SORT_RESULTS = 2 //!< sort points by distances, when k > 1; do not sort by default
		};
		
		//! Find the k nearest neighbours of query
		/*!	\param query query point
		 *	\param k number of nearest neighbour requested
		 *	\param epsilon maximal percentage of error for approximate search, 0 for exact search
		 *	\param optionFlags search options, must be a binary OR of SearchOptionFlags
		 *	\return a vector of k indices to points in cloud */
		virtual IndexVector knn(const Vector& query, const Index k = 1, const T epsilon = 0, const unsigned optionFlags = 0) = 0;
		
		//! Find the k nearest neighbours for each point of query
		/*!	\param query query points
		 *	\param k number of nearest neighbour requested
		 *	\param epsilon maximal percentage of error for approximate search, 0 for exact search
		 *	\param optionFlags search options, must be a binary OR of SearchOptionFlags
		 *	\return a matrix of k x query.cols() indices to points in cloud */
		virtual IndexMatrix knnM(const Matrix& query, const Index k = 1, const T epsilon = 0, const unsigned optionFlags = 0);
		
		//! get a constant version of the statistics
		const Statistics& getStatistics() const { return statistics; }
		
		//! Create a nearest-neighbour search
		/*!	\param cloud data-point cloud in which to search
		 * 	\param preferedType type of search
		 * 	\return an object on which to run nearest neighbour queries */
		static NearestNeighbourSearch* create(const Matrix& cloud, const SearchType preferedType);
		
		//! Create a nearest-neighbour search, using brute-force search, useful for comparison only
		/*!	\param cloud data-point cloud in which to search
		 * 	\return an object on which to run nearest neighbour queries */
		static NearestNeighbourSearch* createBruteForce(const Matrix& cloud);
		
		//! Create a nearest-neighbour search, using a kd-tree with linear heap, good for small k (~up to 30)
		/*!	\param cloud data-point cloud in which to search
		 * 	\return an object on which to run nearest neighbour queries */
		static NearestNeighbourSearch* createKDTreeLinearHeap(const Matrix& cloud);
		
		//! Create a nearest-neighbour search, using a kd-tree with tree heap, good for large k (~from 30)
		/*!	\param cloud data-point cloud in which to search
		 * 	\return an object on which to run nearest neighbour queries */
		static NearestNeighbourSearch* createKDTreeTreeHeap(const Matrix& cloud);
		
	protected:
		//! constructor
		NearestNeighbourSearch(const Matrix& cloud);
		
		//! search statistics
		Statistics statistics;
	};
	
	// Convenience typedefs
	
	//! nearest neighbour search with scalars of type float
	typedef NearestNeighbourSearch<float> NNSearchF;
	//! nearest neighbour search with scalars of type double
	typedef NearestNeighbourSearch<double> NNSearchD;
	
	//@}
}

#endif // __NABO_H
