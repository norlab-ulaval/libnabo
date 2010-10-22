/*
 * 
 C opyr*ight (c) 2010, Stephane Magnenat, ASL, ETHZ, Switzerland
 You can contact the author at <stephane at magnenat dot net>
 
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * Neither the name of the <organization> nor the
 * names of its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * 
 T HIS *SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
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

#ifdef HAVE_OPENCL

#include "CL/cl.hpp"
#include "nabo_private.h"
#include "index_heap.h"
#include <iostream>
#include <stdexcept>
#include <limits>
#include <queue>
#include <algorithm>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>

/*!	\file kdtree_opencl.cpp
 \ *brief kd-tree search, opencl implementation
 \ingroup private
 */

namespace Nabo
{
	//! \ingroup private
	//@{
	
	using namespace std;
	
	template<typename T>
	size_t argMax(const typename NearestNeighbourSearch<T>::Vector& v)
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

	template<typename T>
	size_t KDTreeBalancedPtInLeavesStackOpenCL<T>::getTreeSize(size_t elCount) const
	{
		// FIXME: 64 bits safe stuff, only work for 2^32 elements right now
		assert(elCount > 0);
		elCount --;
		size_t count = 0;
		int i = 31;
		for (; i >= 0; --i)
		{
			if (elCount & (1 << i))
				break;
		}
		for (int j = 0; j <= i; ++j)
			count |= (1 << j);
		count <<= 1;
		count |= 1;
		return count;
	}

	template<typename T>
	void KDTreeBalancedPtInLeavesStackOpenCL<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos, const Vector minValues, const Vector maxValues)
	{
		const size_t count(last - first);
		//cerr << count << endl;
		if (count == 1)
		{
			const int dim = -2-(first->index);
			assert(pos < nodes.size());
			nodes[pos] = Node(dim);
			return;
		}
		
		// find the largest dimension of the box
		size_t cutDim = argMax<T>(maxValues - minValues);
		
		// compute number of elements
		const size_t rightCount(count/2);
		const size_t leftCount(count - rightCount);
		assert(last - rightCount == first + leftCount);
		
		// sort
		nth_element(first, first + leftCount, last, CompareDim(cutDim));
		
		// set node
		const T cutVal((first+leftCount)->pos.coeff(cutDim));
		nodes[pos] = Node(cutDim, cutVal);
		
		//cerr << pos << " cutting on " << cutDim << " at " << (first+leftCount)->pos[cutDim] << endl;
		
		// update bounds for left
		Vector leftMaxValues(maxValues);
		leftMaxValues[cutDim] = cutVal;
		// update bounds for right
		Vector rightMinValues(minValues);
		rightMinValues[cutDim] = cutVal;
		
		// recurse
		buildNodes(first, first + leftCount, childLeft(pos), minValues, leftMaxValues);
		buildNodes(first + leftCount, last, childRight(pos), rightMinValues, maxValues);
	}

	template<typename T>
	KDTreeBalancedPtInLeavesStackOpenCL<T>::KDTreeBalancedPtInLeavesStackOpenCL(const Matrix& cloud, const bool tryGPU):
	NearestNeighbourSearch<T>::NearestNeighbourSearch(cloud)
	{
		// build point vector and compute bounds
		BuildPoints buildPoints;
		buildPoints.reserve(cloud.cols());
		for (int i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			buildPoints.push_back(BuildPoint(v, i));
			const_cast<Vector&>(minBound) = minBound.cwise().min(v);
			const_cast<Vector&>(maxBound) = maxBound.cwise().max(v);
		}
		
		// create nodes
		nodes.resize(getTreeSize(cloud.cols()));
		buildNodes(buildPoints.begin(), buildPoints.end(), 0, minBound, maxBound);
		//for (size_t i = 0; i < nodes.size(); ++i)
		//	cout << i << ": " << nodes[i].dim << " " << nodes[i].cutVal << endl;
	}

	template<typename T>
	typename KDTreeBalancedPtInLeavesStackOpenCL<T>::IndexVector KDTreeBalancedPtInLeavesStackOpenCL<T>::knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T>::ALLOW_SELF_MATCH);
		
		assert(nodes.size() > 0);
		Heap heap(k);
		Vector off(Vector::Zero(query.size()));
		
		statistics.lastQueryVisitCount = 0;
		
		recurseKnn(query, 0, 0, heap, off, 1 + epsilon, allowSelfMatch);
		
		if (optionFlags & NearestNeighbourSearch<T>::SORT_RESULTS)
			heap.sort();
		
		statistics.totalVisitCount += statistics.lastQueryVisitCount;
		
		return heap.getIndexes();
	}

	template<typename T>
	void KDTreeBalancedPtInLeavesStackOpenCL<T>::recurseKnn(const Vector& query, const size_t n, T rd, Heap& heap, Vector& off, const T maxError, const bool allowSelfMatch)
	{
		const Node& node(nodes[n]);
		const int cd(node.dim);
		
		++statistics.lastQueryVisitCount;
		
		if (cd < 0)
		{
			if (cd == -1)
				return;
			const int index(-(cd + 2));
			const T dist(dist2<T>(query, cloud.col(index)));
			if ((dist < heap.headValue()) &&
				(allowSelfMatch || (dist > numeric_limits<T>::epsilon()))
			)
				heap.replaceHead(index, dist);
		}
		else
		{
			const T old_off(off.coeff(cd));
			const T new_off(query.coeff(cd) - node.cutVal);
			if (new_off > 0)
			{
				recurseKnn(query, childRight(n), rd, heap, off, maxError, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd * maxError < heap.headValue())
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childLeft(n), rd, heap, off, maxError, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
			else
			{
				recurseKnn(query, childLeft(n), rd, heap, off, maxError, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd * maxError < heap.headValue())
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childRight(n), rd, heap, off, maxError, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
		}
	}

	template struct KDTreeBalancedPtInLeavesStackOpenCL<float>;
	template struct KDTreeBalancedPtInLeavesStackOpenCL<double>;
	
	//@}
}

#endif // HAVE_OPENCL
