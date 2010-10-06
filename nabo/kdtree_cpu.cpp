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
#include <iostream>
#include <stdexcept>
#include <limits>
#include <queue>
#include <algorithm>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>

namespace Nabo
{
	using namespace std;
	
	template<typename T>
	size_t argMax(const typename NearestNeighborSearch<T>::Vector& v)
	{
		T maxVal(0);
		size_t maxIdx(0);
		for (int i = 0; i < v.size(); ++i)
		{
			if (v[i] > maxVal)
			{
				maxVal = v[i];
				maxIdx = i;
			}
		}
		return maxIdx;
	}
	
	// OPT
	template<typename T, typename Heap>
	pair<T,T> KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::getBounds(const BuildPointsIt first, const BuildPointsIt last, const unsigned dim)
	{
		T minVal(boost::numeric::bounds<T>::highest());
		T maxVal(boost::numeric::bounds<T>::lowest());
		
		for (BuildPointsCstIt it(first); it != last; ++it)
		{
			const T val(cloud.coeff(dim, *it));
			minVal = min(val, minVal);
			maxVal = max(val, maxVal);
		}
		
		return make_pair<T>(minVal, maxVal);
	}
	
	template<typename T, typename Heap>
	unsigned KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const Vector minValues, const Vector maxValues)
	{
		const int count(last - first);
		const unsigned pos(nodes.size());
		
		//cerr << count << endl;
		if (count == 1)
		{
			nodes.push_back(Node(*first, &cloud.coeff(0, *first)));
			return pos;
		}
		
		// find the largest dimension of the box
		const unsigned cutDim = argMax<T>(maxValues - minValues);
		const T idealCutVal((maxValues(cutDim) + minValues(cutDim))/2);
		
		// get bounds from actual points
		const pair<T,T> minMaxVals(getBounds(first, last, cutDim));
		
		// correct cut following bounds
		T cutVal;
		if (idealCutVal < minMaxVals.first)
			cutVal = minMaxVals.first;
		else if (idealCutVal > minMaxVals.second)
			cutVal = minMaxVals.second;
		else
			cutVal = idealCutVal;
		
		int l(0);
		int r(count-1);
		// partition points around cutVal
		while (1)
		{
			while (l < count && cloud.coeff(cutDim, *(first+l)) < cutVal)
				++l;
			while (r >= 0 && cloud.coeff(cutDim, *(first+r)) >= cutVal)
				--r;
			if (l > r)
				break;
			swap(*(first+l), *(first+r));
			++l; --r;
		}
		const int br1 = l;	// now: points[0..br1-1] < cutVal <= points[br1..count-1]
		r = count-1;
		// partition points[br1..count-1] around cv
		while (1)
		{
			while (l < count && cloud.coeff(cutDim, *(first+l)) <= cutVal)
				++l;
			while (r >= br1 && cloud.coeff(cutDim, *(first+r)) > cutVal)
				--r;
			if (l > r)
				break;
			swap(*(first+l), *(first+r));
			++l; --r;
		}
		const int br2 = l; // now: points[br1..br2-1] == cv < points[br2..count-1]
		
		// find best split index
		int leftCount;
		if (idealCutVal < minMaxVals.first)
			leftCount = 1;
		else if (idealCutVal > minMaxVals.second)
			leftCount = count-1;
		else if (br1 > count / 2)
			leftCount = br1;
		else if (br2 < count / 2)
			leftCount = br2;
		else
			leftCount = count / 2;
		assert(leftCount > 0);
		assert(leftCount < count);
		
		// update bounds for left
		Vector leftMaxValues(maxValues);
		leftMaxValues[cutDim] = cutVal;
		// update bounds for right
		Vector rightMinValues(minValues);
		rightMinValues[cutDim] = cutVal;
		
		// add this
		nodes.push_back(Node(cutDim, cutVal, 0));
		
		// recurse
		const unsigned __attribute__ ((unused)) leftChild = buildNodes(first, first + leftCount, minValues, leftMaxValues);
		assert(leftChild == pos + 1);
		const unsigned rightChild = buildNodes(first + leftCount, last, rightMinValues, maxValues);
		
		// write right child index and return
		nodes[pos].rightChild = rightChild;
		return pos;
	}

	template<typename T, typename Heap>
	KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt(const Matrix& cloud):
		NearestNeighborSearch<T>::NearestNeighborSearch(cloud),
		dimCount(cloud.rows())
	{
		// build point vector and compute bounds
		BuildPoints buildPoints;
		buildPoints.reserve(cloud.cols());
		for (int i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			buildPoints.push_back(i);
			const_cast<Vector&>(minBound) = minBound.cwise().min(v);
			const_cast<Vector&>(maxBound) = maxBound.cwise().max(v);
		}
		
		// create nodes
		//nodes.resize(getTreeSize(cloud.cols()));
		buildNodes(buildPoints.begin(), buildPoints.end(), minBound, maxBound);
		buildPoints.clear();
		//for (size_t i = 0; i < nodes.size(); ++i)
		//	cout << i << ": " << nodes[i].dim << " " << nodes[i].cutVal <<  " " << nodes[i].rightChild << endl;
		
	}
	
	template<typename T, typename Heap>
	typename KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::IndexVector KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighborSearch<T>::ALLOW_SELF_MATCH);
		
		assert(nodes.size() > 0);
		Heap heap(k);
		std::vector<T> off(dimCount, 0);
		
		statistics.lastQueryVisitCount = 0;
		
		if (allowSelfMatch)
			recurseKnn<true>(&query.coeff(0), 0, 0, heap, off, 1+epsilon);
		else
			recurseKnn<false>(&query.coeff(0), 0, 0, heap, off, 1+epsilon);
		
		if (optionFlags & NearestNeighborSearch<T>::SORT_RESULTS)
			heap.sort();
		
		statistics.totalVisitCount += statistics.lastQueryVisitCount;
		
		return heap.getIndexes();
	}
	
	
	template<typename T, typename Heap>
	typename KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::IndexMatrix KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::knnM(const Matrix& query, const Index k, const T epsilon, const unsigned optionFlags) 
	{
		const bool allowSelfMatch(optionFlags & NearestNeighborSearch<T>::ALLOW_SELF_MATCH);
		assert(nodes.size() > 0);
		
		assert(nodes.size() > 0);
		Heap heap(k);
		
		std::vector<T> off(dimCount, 0);
		
		IndexMatrix result(k, query.cols());
		const int colCount(query.cols());
		for (int i = 0; i < colCount; ++i)
		{
			fill(off.begin(), off.end(), 0);
			heap.reset();
			
			// FIXME: add define for statistics
			statistics.lastQueryVisitCount = 0;
			
			if (allowSelfMatch)
				recurseKnn<true>(&query.coeff(0, i), 0, 0, heap, off, 1+epsilon);
			else
				recurseKnn<false>(&query.coeff(0, i), 0, 0, heap, off, 1+epsilon);
			
			if (optionFlags & NearestNeighborSearch<T>::SORT_RESULTS)
				heap.sort();
			
			result.col(i) = heap.getIndexes();
			
			statistics.totalVisitCount += statistics.lastQueryVisitCount;
		}
		
		return result;
	}
	
	template<typename T, typename Heap> template<bool allowSelfMatch>
	void KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap>::recurseKnn(const T* query, const unsigned n, T rd, Heap& heap, std::vector<T>& off, const T maxError)
	{
		const Node& node(nodes[n]);
		//++statistics.lastQueryVisitCount;
		
		if (node.rightChild == Node::INVALID_CHILD)
		{
			//const T dist(dist2<T>(query, cloud.col(index)));
			//const T dist((query - cloud.col(index)).squaredNorm());
			T dist(0);
			const T* qPtr(query);
			const T* dPtr(node.pt);
			for (int i = 0; i < dimCount; ++i)
			{
				const T diff(*qPtr - *dPtr);
				dist += diff*diff;
				qPtr++; dPtr++;
			}
			if ((dist < heap.headValue()) &&
				(allowSelfMatch || (dist > numeric_limits<T>::epsilon()))
			)
				heap.replaceHead(node.dim, dist);
		}
		else
		{
			const Index cd(node.dim);
			T& offcd(off[cd]);
			//const T old_off(off.coeff(cd));
			const T old_off(offcd);
			const T new_off(query[cd] - node.cutVal);
			if (new_off > 0)
			{
				recurseKnn<allowSelfMatch>(query, node.rightChild, rd, heap, off, maxError);
				rd += - old_off*old_off + new_off*new_off;
				if (rd * maxError < heap.headValue())
				{
					offcd = new_off;
					recurseKnn<allowSelfMatch>(query, n + 1, rd, heap, off, maxError);
					offcd = old_off;
				}
			}
			else
			{
				recurseKnn<allowSelfMatch>(query, n+1, rd, heap, off, maxError);
				rd += - old_off*old_off + new_off*new_off;
				if (rd * maxError < heap.headValue())
				{
					offcd = new_off;
					recurseKnn<allowSelfMatch>(query, node.rightChild, rd, heap, off, maxError);
					offcd = old_off;
				}
			}
		}
	}
	
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapSTL<int,float>>;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapBruteForceVector<int,float>>;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapSTL<int,double>>;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapBruteForceVector<int,double>>;
}
