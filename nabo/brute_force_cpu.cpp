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

#include "nabo_private.h"
#include "index_heap.h"

/*!	\file brute_force_cpu.cpp
	\brief brute force search, cpu implementation
	\ingroup private
*/

namespace Nabo
{
	using namespace std;
	
	template<typename T>
	BruteForceSearch<T>::BruteForceSearch(const Matrix& cloud):
		NearestNeighbourSearch<T>::NearestNeighbourSearch(cloud)
	{
		// compute bounds
		for (int i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			const_cast<Vector&>(this->minBound) = this->minBound.cwise().min(v);
			const_cast<Vector&>(this->maxBound) = this->maxBound.cwise().max(v);
		}
	}
	

	template<typename T>
	typename NearestNeighbourSearch<T>::IndexVector BruteForceSearch<T>::knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T>::ALLOW_SELF_MATCH);
		
		IndexHeapSTL<Index, T> heap(k);
		
		for (int i = 0; i < this->cloud.cols(); ++i)
		{
			const T dist(dist2<T>(this->cloud.col(i), query));
			if ((dist < heap.headValue()) &&
				(allowSelfMatch || (dist > numeric_limits<T>::epsilon())))
				heap.replaceHead(i, dist);
		}
		
		this->statistics.lastQueryVisitCount = this->cloud.cols();
		this->statistics.totalVisitCount += this->statistics.lastQueryVisitCount;
		
		if (optionFlags & NearestNeighbourSearch<T>::SORT_RESULTS)
			heap.sort();
		
		return heap.getIndexes();
	}
	
	template struct BruteForceSearch<float>;
	template struct BruteForceSearch<double>;
}
