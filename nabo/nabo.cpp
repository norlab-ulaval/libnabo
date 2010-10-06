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

namespace Nabo
{
	using namespace std;
	
	template<typename T>
	NearestNeighborSearch<T>::NearestNeighborSearch(const Matrix& cloud):
		cloud(cloud),
		dim(cloud.rows()),
		minBound(Vector::Constant(dim, numeric_limits<T>::max())),
		maxBound(Vector::Constant(dim, numeric_limits<T>::min()))
	{
		
	}
	
	template<typename T>
	typename NearestNeighborSearch<T>::IndexMatrix NearestNeighborSearch<T>::knnM(const Matrix& query, const Index k, const T epsilon, const unsigned optionFlags) 
	{
		IndexMatrix result(k, query.cols());
		const int colCount(query.cols());
		
		for (int i = 0; i < colCount; ++i)
		{
			const Vector& q(query.col(i));
			result.col(i) = knn(q, k, epsilon, optionFlags);
		}
		
		return result;
	}
	
	template<typename T>
	NearestNeighborSearch<T>* NearestNeighborSearch<T>::create(const Matrix& cloud, const SearchType preferedType)
	{
		switch (preferedType)
		{
			case BRUTE_FORCE: return new BruteForceSearch<T>(cloud);
			case KDTREE_LINEAR_HEAP: return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapBruteForceVector<int,T>>(cloud);
			case KDTREE_TREE_HEAP: return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapSTL<int,T>>(cloud);
			default: throw runtime_error("Unknown search type");
		}
	}
	
	template<typename T>
	NearestNeighborSearch<T>* NearestNeighborSearch<T>::createBruteForce(const Matrix& cloud)
	{
		return new BruteForceSearch<T>(cloud);;
	}
	
	template<typename T>
	NearestNeighborSearch<T>* NearestNeighborSearch<T>::createKDTreeLinearHeap(const Matrix& cloud)
	{
		return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapBruteForceVector<int,T>>(cloud);
	}
	
	template<typename T>
	NearestNeighborSearch<T>* NearestNeighborSearch<T>::createKDTreeTreeHeap(const Matrix& cloud)
	{
		return new KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, IndexHeapSTL<int,T>>(cloud);
	}
	
	template struct NearestNeighborSearch<float>;
	template struct NearestNeighborSearch<double>;
}
