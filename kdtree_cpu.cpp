#include "nabo.h"
#include "index_heap.h"
#include <iostream>
#include <stdexcept>
#include <limits>
#include <queue>
#include <algorithm>

namespace Nabo
{
	using namespace std;

	template<typename T>
	size_t KDTree<T>::getTreeSize(size_t elCount) const
	{
		// FIXME: 64 bits safe stuff, only work for 2^32 elements right now
		size_t count = 0;
		int i = 31;
		for (; i >= 0; --i)
		{
			if (elCount & (1 << i))
				break;
		}
		for (int j = 0; j <= i; ++j)
			count |= (1 << j);
		//cerr << "tree size " << count << " (" << elCount << " elements)\n";
		return count;
	}

	template<typename T>
	size_t KDTree<T>::argMax(const Vector& v) const
	{
		T maxVal(0);
		size_t maxIdx;
		for (size_t i = 0; i < this->dim; ++i)
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
	typename KDTree<T>::Indexes KDTree<T>::cloudIndexesFromNodesIndexes(const Indexes& indexes) const
	{
		Indexes cloudIndexes;
		cloudIndexes.reserve(indexes.size());
		for (typename Indexes::const_iterator it(indexes.begin()); it != indexes.end(); ++it)
			cloudIndexes.push_back(nodes[*it].index);
		return cloudIndexes;
	}

	template<typename T>
	void KDTree<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos)
	{
		const size_t count(last - first);
		//cerr << count << endl;
		if (count == 1)
		{
			nodes[pos] = Node(first->pos, -1, first->index);
			return;
		}
		
		// estimate variance
		// get mean
		Vector mean(Vector::Zero(this->dim));
		for (BuildPointsCstIt it(first); it != last; ++it)
			mean += it->pos;
		mean /= last - first;
		// get sum of variance
		Vector var(Vector::Zero(this->dim));
		for (BuildPointsCstIt it(first); it != last; ++it)
			var += (it->pos - mean).cwise() * (it->pos - mean);
		// get dimension of maxmial variance
		const size_t cutDim = argMax(var);
		
		// sort
		sort(first, last, CompareDim(cutDim));
		
		// set node
		const size_t recurseCount(count-1);
		const size_t rightCount(recurseCount/2);
		const size_t leftCount(recurseCount-rightCount);
		assert(last - rightCount == first + leftCount + 1);
		
		nodes[pos] = Node((first+leftCount)->pos, cutDim, (first+leftCount)->index);
		
		//cerr << pos << " cutting on " << cutDim << " at " << (first+leftCount)->pos[cutDim] << endl;
		
		// recurse
		if (count > 2)
		{
			buildNodes(first, first + leftCount, childLeft(pos));
			buildNodes(first + leftCount + 1, last, childRight(pos));
		}
		else
		{
			nodes[childLeft(pos)] = Node(first->pos, -1, first->index);
			nodes[childRight(pos)] = Node(Vector(), -2, 0);
		}
	}

	template<typename T>
	void KDTree<T>::dump(const Vector minValues, const Vector maxValues, const size_t pos) const
	{
		const Node& node(nodes[pos]);
		
		if (node.dim >= -1)
		{
			if (this->dim == 2)
				cout << "<circle cx=\"" << 100*node.pos(0) << "\" cy=\"" << 100*node.pos(1) << "\" r=\"1\" stroke=\"black\" stroke-width=\"0.2\" fill=\"red\"/>" << endl;
			else
				cout << "pt at\n" << node.pos << endl;
		}
		if (node.dim >= 0)
		{
			//cerr << "in bounds:\n" << minValues << "\nto\n" << maxValues << endl;
			
			// update bounds for left
			Vector leftMaxValues(maxValues);
			leftMaxValues[node.dim] = node.pos[node.dim];
			// update bounds for right
			Vector rightMinValues(minValues);
			rightMinValues[node.dim] = node.pos[node.dim];
			
			// print line
			if (this->dim == 2)
				cout << "<line x1=\"" << 100*rightMinValues(0) << "\" y1=\"" << 100*rightMinValues(1) << "\" x2=\"" << 100*leftMaxValues(0) << "\" y2=\"" << 100*leftMaxValues(1) << "\" style=\"stroke:rgb(0,0,0);stroke-width:0.2\"/>" << endl;
			else
				cout << "cut from\n" << rightMinValues << "\nto\n" << leftMaxValues << endl;
			// recurs
			dump(minValues, leftMaxValues, childLeft(pos));
			dump(rightMinValues, maxValues, childRight(pos));
		}
	}

	template<typename T>
	KDTree<T>::KDTree(const Matrix& cloud):
		NearestNeighborSearch<T>::NearestNeighborSearch(cloud)
	{
		// build point vector and compute bounds
		vector<BuildPoint> buildPoints;
		for (size_t i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			buildPoints.push_back(BuildPoint(v, i));
			const_cast<Vector&>(this->minBound) = this->minBound.cwise().min(v);
			const_cast<Vector&>(this->maxBound) = this->maxBound.cwise().max(v);
		}
		
		// create nodes
		nodes.resize(getTreeSize(cloud.cols()));
		buildNodes(buildPoints.begin(), buildPoints.end(), 0);
		
		// dump nodes
		//dump(minBound, maxBound, 0);
	}

	static int totKDTreeVisitCount = 0;

	template<typename T>
	typename KDTree<T>::Indexes KDTree<T>::knn(const Vector& query, const Index k, const bool allowSelfMatch)
	{
		typedef priority_queue<SearchElement> Queue;
		
		Queue queue;
		queue.push(SearchElement(0, 0));
		IndexHeap<Index, T> heap(k);
		this->statistics.lastQueryVisitCount = 0;
		
		while (!queue.empty())
		{
			SearchElement el(queue.top());
			queue.pop();
			
			// nothing is closer, we found best
			if (el.minDist > heap.head().value)
				break;
			
			size_t n(el.index);
			while (1)
			{
				const Node& node(nodes[n]);
				assert (node.dim != -2);
				
				const T dist(dist2(node.pos, query));
				if ((dist < heap.head().value) &&
					(allowSelfMatch || (dist > numeric_limits<T>::epsilon())))
					heap.replaceHead(n, dist);
				
				// if we are at leaf, stop
				if (node.dim < 0)
					break;
				
				const T offset(query[node.dim] - node.pos[node.dim]);
				const T offset2(offset * offset);
				const T bestDist(heap.head().value);
				if (offset > 0)
				{
					// enqueue offside ?
					if (offset2 < bestDist && nodes[childLeft(n)].dim != -2)
						queue.push(SearchElement(childLeft(n), offset2));
					// continue onside
					if (nodes[childRight(n)].dim != -2)
						n = childRight(n);
					else
						break;
				}
				else
				{
					// enqueue offside ?
					if (offset2 < bestDist && nodes[childRight(n)].dim != -2)
						queue.push(SearchElement(childRight(n), offset2));
					// continue onside
					if (nodes[childLeft(n)].dim != -2)
						n = childLeft(n);
					else
						break;
				}
				++this->statistics.lastQueryVisitCount;
			}
		}
		this->statistics.totalVisitCount += this->statistics.lastQueryVisitCount;
		
		// TODO: add flag
		heap.sort();
		
		return cloudIndexesFromNodesIndexes(heap.getIndexes());
	}
	
	template struct KDTree<float>;
	template struct KDTree<double>;
}
