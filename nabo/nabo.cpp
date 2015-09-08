/*

Copyright (c) 2010--2011, Stephane Magnenat, ASL, ETHZ, Switzerland
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

	template<typename T, typename CloudType>
	NearestNeighbourSearch<T, CloudType>::NearestNeighbourSearch(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags):
		cloud(cloud),
		dim(min(dim, int(cloud.rows()))),
		creationOptionFlags(creationOptionFlags),
		minBound(Vector::Constant(this->dim, numeric_limits<T>::max())),
		maxBound(Vector::Constant(this->dim, numeric_limits<T>::min()))
	{
		if (cloud.cols() == 0)
			throw runtime_error("Cloud has no points");
		if (cloud.rows() == 0)
			throw runtime_error("Cloud has 0 dimensions");
	}
	
	template<typename T, typename CloudType>
	unsigned long NearestNeighbourSearch<T, CloudType>::knn(const Vector& query, IndexVector& indices, Vector& dists2, const Index k, const T epsilon, const unsigned optionFlags, const T maxRadius) const
	{
#ifdef EIGEN3_API
		const Eigen::Map<const Matrix> queryMatrix(&query.coeff(0,0), dim, 1);
#else // EIGEN3_API
		const Eigen::Map<Matrix> queryMatrix(&query.coeff(0,0), dim, 1);
#endif // EIGEN3_API
		// note: this is inefficient, because we copy memory, due to the template-
		// based abstraction of Eigen. High-performance implementation should
		// take care of knnM and then implement knn on top of it.
		// C++0x should solve this with rvalue
		IndexMatrix indexMatrix(k, 1);
		Matrix dists2Matrix(k, 1);
		const unsigned long stats = knn(queryMatrix, indexMatrix, dists2Matrix, k, epsilon, optionFlags, maxRadius);
		indices = indexMatrix.col(0);
		dists2 = dists2Matrix.col(0);
		return stats;
	}

	template<typename T, typename CloudType>
	void NearestNeighbourSearch<T, CloudType>::checkSizesKnn(const Matrix& query, const IndexMatrix& indices, const Matrix& dists2, const Index k, const unsigned optionFlags, const Vector* maxRadii) const
	{
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T, CloudType>::ALLOW_SELF_MATCH);
		if (allowSelfMatch)
		{
			if (k > cloud.cols())
				throw runtime_error((boost::format("Requesting more points (%1%) than available in cloud (%2%)") % k % cloud.cols()).str());
		}
		else
		{
			if (k > cloud.cols()-1)
				throw runtime_error((boost::format("Requesting more points (%1%) than available in cloud minus 1 (%2%) (as self match is forbidden)") % k % (cloud.cols()-1)).str());
		}
		if (query.rows() < dim)
			throw runtime_error((boost::format("Query has less dimensions (%1%) than requested for cloud (%2%)") % query.rows() % dim).str());
		if (indices.rows() != k)
			throw runtime_error((boost::format("Index matrix has a different number of rows (%1%) than k (%2%)") % indices.rows() % k).str());
		if (indices.cols() != query.cols())
			throw runtime_error((boost::format("Index matrix has a different number of columns (%1%) than query (%2%)") % indices.rows() % query.cols()).str());
		if (dists2.rows() != k)
			throw runtime_error((boost::format("Distance matrix has a different number of rows (%1%) than k (%2%)") % dists2.rows() % k).str());
		if (dists2.cols() != query.cols())
			throw runtime_error((boost::format("Distance matrix has a different number of columns (%1%) than query (%2%)") % dists2.rows() % query.cols()).str());
		if (maxRadii && (maxRadii->size() != query.cols()))
			throw runtime_error((boost::format("Maximum radii vector has not the same length (%1%) than query has columns (%2%)") % maxRadii->size() % k).str());
		const unsigned maxOptionFlagsValue(ALLOW_SELF_MATCH|SORT_RESULTS);
		if (optionFlags > maxOptionFlagsValue)
			throw runtime_error((boost::format("OR-ed value of option flags (%1%) is larger than maximal valid value (%2%)") % optionFlags % maxOptionFlagsValue).str());
	}
	

	template<typename T, typename CloudType>
	NearestNeighbourSearch<T, CloudType>* NearestNeighbourSearch<T, CloudType>::create(const CloudType& cloud, const Index dim, const SearchType preferedType, const unsigned creationOptionFlags, const Parameters& additionalParameters)
	{
		if (dim <= 0)
			throw runtime_error("Your space must have at least one dimension");
		switch (preferedType)
		{
			case BRUTE_FORCE: return new BruteForceSearch<T, CloudType>(cloud, dim, creationOptionFlags);
			case KDTREE_LINEAR_HEAP: return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapBruteForceVector<int,T>, CloudType>(cloud, dim, creationOptionFlags, additionalParameters);
			case KDTREE_TREE_HEAP: return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapSTL<int,T>, CloudType>(cloud, dim, creationOptionFlags, additionalParameters);
			#ifdef HAVE_OPENCL
			case KDTREE_CL_PT_IN_NODES: return new KDTreeBalancedPtInNodesStackOpenCL<T, CloudType>(cloud, dim, creationOptionFlags, CL_DEVICE_TYPE_GPU);
			case KDTREE_CL_PT_IN_LEAVES: return new KDTreeBalancedPtInLeavesStackOpenCL<T, CloudType>(cloud, dim, creationOptionFlags, CL_DEVICE_TYPE_GPU);
			case BRUTE_FORCE_CL: return new BruteForceSearchOpenCL<T, CloudType>(cloud, dim, creationOptionFlags, CL_DEVICE_TYPE_GPU);
			#else // HAVE_OPENCL
			case KDTREE_CL_PT_IN_NODES: throw runtime_error("OpenCL not found during compilation");
			case KDTREE_CL_PT_IN_LEAVES: throw runtime_error("OpenCL not found during compilation");
			case BRUTE_FORCE_CL: throw runtime_error("OpenCL not found during compilation");
			#endif // HAVE_OPENCL
			default: throw runtime_error("Unknown search type");
		}
	}

	template<typename T, typename CloudType>
	NearestNeighbourSearch<T, CloudType>* NearestNeighbourSearch<T, CloudType>::createBruteForce(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags)
	{
		if (dim <= 0)
			throw runtime_error("Your space must have at least one dimension");
		return new BruteForceSearch<T, CloudType>(cloud, dim, creationOptionFlags);
	}

	template<typename T, typename CloudType>
	NearestNeighbourSearch<T, CloudType>* NearestNeighbourSearch<T, CloudType>::createKDTreeLinearHeap(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags, const Parameters& additionalParameters)
	{
		if (dim <= 0)
			throw runtime_error("Your space must have at least one dimension");
		return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapBruteForceVector<int,T>, CloudType>(cloud, dim, creationOptionFlags, additionalParameters);
	}

	template<typename T, typename CloudType>
	NearestNeighbourSearch<T, CloudType>* NearestNeighbourSearch<T, CloudType>::createKDTreeTreeHeap(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags, const Parameters& additionalParameters)
	{
		if (dim <= 0)
			throw runtime_error("Your space must have at least one dimension");
		return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapSTL<int,T>, CloudType>(cloud, dim, creationOptionFlags, additionalParameters);
	}
	
	template struct NearestNeighbourSearch<float>;
	template struct NearestNeighbourSearch<double>;
	template struct NearestNeighbourSearch<float, Eigen::Matrix3Xf>;
	template struct NearestNeighbourSearch<double, Eigen::Matrix3Xd>;
	template struct NearestNeighbourSearch<float, Eigen::Map<const Eigen::Matrix3Xf, Eigen::Aligned> >;
	template struct NearestNeighbourSearch<double, Eigen::Map<const Eigen::Matrix3Xd, Eigen::Aligned> >;
}
