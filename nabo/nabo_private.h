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

#ifndef __NABO_PRIVATE_H
#define __NABO_PRIVATE_H

#include "nabo.h"

namespace Nabo
{
	// Euclidean distance
	template<typename T, typename A, typename B>
	inline T dist2(const A& v0, const B& v1)
	{
		return (v0 - v1).squaredNorm();
	}

	// Brute-force nearest neighbor
	template<typename T>
	struct BruteForceSearch: public NearestNeighborSearch<T>
	{
		typedef typename NearestNeighborSearch<T>::Vector Vector;
		typedef typename NearestNeighborSearch<T>::Matrix Matrix;
		typedef typename NearestNeighborSearch<T>::Index Index;
		typedef typename NearestNeighborSearch<T>::IndexVector IndexVector;

		BruteForceSearch(const Matrix& cloud);
		virtual IndexVector knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags);
	};
	
	//  KDTree, unbalanced, points in leaves, stack, implicit bounds, ANN_KD_SL_MIDPT, optimised
	template<typename T, typename Heap>
	struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt: public NearestNeighborSearch<T>
	{
		typedef typename NearestNeighborSearch<T>::Vector Vector;
		typedef typename NearestNeighborSearch<T>::Matrix Matrix;
		typedef typename NearestNeighborSearch<T>::Index Index;
		typedef typename NearestNeighborSearch<T>::IndexVector IndexVector;
		typedef typename NearestNeighborSearch<T>::IndexMatrix IndexMatrix;
		
		using NearestNeighborSearch<T>::statistics;
		using NearestNeighborSearch<T>::cloud;
		using NearestNeighborSearch<T>::minBound;
		using NearestNeighborSearch<T>::maxBound;
		
	protected:
		typedef std::vector<Index> BuildPoints;
		typedef typename BuildPoints::iterator BuildPointsIt;
		typedef typename BuildPoints::const_iterator BuildPointsCstIt;
		
		struct Node
		{
			enum
			{
				INVALID_CHILD = 0xffffffff,
				INVALID_PT = 0
			};
			Index dim; // also index for point
			unsigned rightChild;
			union
			{
				T cutVal;
				const T* pt;
			};
			
			Node(const Index dim, const T cutVal, unsigned rightChild):
				dim(dim), rightChild(rightChild), cutVal(cutVal) {}
			Node(const Index index = 0, const T* pt = 0):
				dim(index), rightChild(INVALID_CHILD), pt(pt) {}
		};
		typedef std::vector<Node> Nodes;
		
		Nodes nodes;
		const int dimCount;
		
		std::pair<T,T> getBounds(const BuildPointsIt first, const BuildPointsIt last, const unsigned dim);
		unsigned buildNodes(const BuildPointsIt first, const BuildPointsIt last, const Vector minValues, const Vector maxValues);
		
		template<bool allowSelfMatch>
		void recurseKnn(const T* query, const unsigned n, T rd, Heap& heap, std::vector<T>& off, const T maxError);
		
	public:
		KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt(const Matrix& cloud);
		virtual IndexVector knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags);
		virtual IndexMatrix knnM(const Matrix& query, const Index k, const T epsilon, const unsigned optionFlags);
	};
}

#endif // __NABO_H
