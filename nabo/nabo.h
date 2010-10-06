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

namespace Nabo
{
	// Nearest neighbor search interface, templatized on scalar type
	template<typename T>
	struct NearestNeighborSearch
	{
		typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
		typedef typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix; // each entry is a col, matrix has dim rows
		typedef int Index;
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> IndexMatrix;
		
		const Matrix& cloud;
		const size_t dim;
		const Vector minBound;
		const Vector maxBound;
		
		struct Statistics
		{
			Statistics():lastQueryVisitCount(0),totalVisitCount(0) {}
			std::atomic_uint lastQueryVisitCount;
			std::atomic_uint totalVisitCount;
		};
		
		enum SearchType
		{
			BRUTE_FORCE = 0,
			KDTREE_LINEAR_HEAP = 1,
			KDTREE_TREE_HEAP = 2,
			SEARCH_TYPE_COUNT
		};
		
		enum SearchOptionFlags
		{
			ALLOW_SELF_MATCH = 1,
			SORT_RESULTS = 2
		};
		
		NearestNeighborSearch(const Matrix& cloud);
		virtual IndexVector knn(const Vector& query, const Index k = 1, const T epsilon = 0, const unsigned optionFlags = 0) = 0;
		virtual IndexMatrix knnM(const Matrix& query, const Index k = 1, const T epsilon = 0, const unsigned optionFlags = 0);
		const Statistics& getStatistics() const { return statistics; }
		
		static NearestNeighborSearch* create(const Matrix& cloud, const SearchType preferedType);
		static NearestNeighborSearch* createBruteForce(const Matrix& cloud);
		static NearestNeighborSearch* createKDTreeLinearHeap(const Matrix& cloud);
		static NearestNeighborSearch* createKDTreeTreeHeap(const Matrix& cloud);
		
	protected:
		Statistics statistics;
	};
	
	// Convenience typedefs
	typedef NearestNeighborSearch<float> NNSearchF;
	typedef NearestNeighborSearch<double> NNSearchD;
}

#endif // __NABO_H
