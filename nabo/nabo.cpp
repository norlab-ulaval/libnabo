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

#include "nabo.h"
#include "nabo_private.h"
#include "index_heap.h"
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <boost/format.hpp>

/*!	\file nabo.cpp
	\brief implementation of public interface
	\ingroup private
*/

namespace Nabo
{
	using namespace std;
	
	template<typename T>
	NearestNeighbourSearch<T>::NearestNeighbourSearch(const Matrix& cloud, const Index dim):
		cloud(cloud),
		dim(min(dim, cloud.rows())),
		minBound(Vector::Constant(this->dim, numeric_limits<T>::max())),
		maxBound(Vector::Constant(this->dim, numeric_limits<T>::min()))
	{
		
	}
	
	template<typename T>
	void NearestNeighbourSearch<T>::knn(const Vector& query, IndexVector& indices, Vector& dists2, const Index k, const T epsilon, const unsigned optionFlags)
	{
		const Eigen::Map<Matrix> queryMatrix(&query.coeff(0,0), dim, 1);
		// note: this is inefficient, because we copy memory, due to the template-
		// based abstraction of Eigen. High-performance implementation should
		// take care of knnM and then implement knn on top of it.
		// C++0x should solve this with rvalue
		IndexMatrix indexMatrix(k, 1);
		Matrix dists2Matrix(k, 1);
		knn(queryMatrix, indexMatrix, dists2Matrix, k, epsilon, optionFlags);
		indices = indexMatrix.col(0);
		dists2 = dists2Matrix.col(0);
	}
	
	template<typename T>
	void NearestNeighbourSearch<T>::checkSizesKnn(const Matrix& query, const IndexMatrix& indices, const Matrix& dists2, const Index k)
	{
		if (query.rows() < dim)
			throw runtime_error((boost::format("Query has less dimensions (%1%) than requested for cloud (%2%)") % query.rows() % dim).str());
		if (indices.rows() != k)
			throw runtime_error((boost::format("Index matrix has less rows (%1%) than k (%2%)") % indices.rows() % k).str());
		if (indices.cols() != query.cols())
			throw runtime_error((boost::format("Index matrix has less columns (%1%) than query (%2%)") % indices.rows() % query.cols()).str());
		if (dists2.rows() != k)
			throw runtime_error((boost::format("Distance matrix has less rows (%1%) than k (%2%)") % dists2.rows() % k).str());
		if (dists2.cols() != query.cols())
			throw runtime_error((boost::format("Distance matrix has less columns (%1%) than query (%2%)") % dists2.rows() % query.cols()).str());
	}
	
	
	template<typename T>
	NearestNeighbourSearch<T>* NearestNeighbourSearch<T>::create(const Matrix& cloud, const Index dim, const SearchType preferedType)
	{
		switch (preferedType)
		{
			case BRUTE_FORCE: return new BruteForceSearch<T>(cloud, dim);
			case KDTREE_LINEAR_HEAP: return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapBruteForceVector<int,T>>(cloud, dim);
			case KDTREE_TREE_HEAP: return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapSTL<int,T>>(cloud, dim);
			#ifdef HAVE_OPENCL
			case KDTREE_CL: return new KDTreeBalancedPtInLeavesStackOpenCL<T>(cloud, dim, CL_DEVICE_TYPE_GPU);
			case BRUTE_FORCE_CL: return new BruteForceSearchOpenCL<T>(cloud, dim, CL_DEVICE_TYPE_GPU);
			#else // HAVE_OPENCL
			case KDTREE_CL_CPU: throw runtime_error("OpenCL not found during compilation");
			case KDTREE_CL_GPU: throw runtime_error("OpenCL not found during compilation");
			case BRUTE_FORCE_CL_CPU: throw runtime_error("OpenCL not found during compilation");
			case BRUTE_FORCE_CL_GPU: throw runtime_error("OpenCL not found during compilation");
			#endif // HAVE_OPENCL
			default: throw runtime_error("Unknown search type");
		}
	}
	
	template<typename T>
	NearestNeighbourSearch<T>* NearestNeighbourSearch<T>::createBruteForce(const Matrix& cloud, const Index dim)
	{
		return new BruteForceSearch<T>(cloud, dim);
	}
	
	template<typename T>
	NearestNeighbourSearch<T>* NearestNeighbourSearch<T>::createKDTreeLinearHeap(const Matrix& cloud, const Index dim)
	{
		return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapBruteForceVector<int,T>>(cloud, dim);
	}
	
	template<typename T>
	NearestNeighbourSearch<T>* NearestNeighbourSearch<T>::createKDTreeTreeHeap(const Matrix& cloud, const Index dim)
	{
		return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapSTL<int,T>>(cloud, dim);
	}
	
	template struct NearestNeighbourSearch<float>;
	template struct NearestNeighbourSearch<double>;
}
