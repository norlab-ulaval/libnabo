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
	size_t argMax(const typename NearestNeighborSearch<T>::Vector& v)
	{
		T maxVal(0);
		size_t maxIdx(0);
		for (size_t i = 0; i < v.size(); ++i)
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
	size_t KDTreePtInNodes<T>::getTreeSize(size_t elCount) const
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
	typename KDTreePtInNodes<T>::IndexVector KDTreePtInNodes<T>::cloudIndexesFromNodesIndexes(const IndexVector& indexes) const
	{
		IndexVector cloudIndexes(indexes.size());
		for (size_t i = 0; i < indexes.size(); ++i)
			cloudIndexes.coeffRef(i) = nodes[indexes[i]].index;
		return cloudIndexes;
	}

	template<typename T>
	void KDTreePtInNodes<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos)
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
		const size_t cutDim = argMax<T>(var);
		
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
	void KDTreePtInNodes<T>::dump(const Vector minValues, const Vector maxValues, const size_t pos) const
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
	KDTreePtInNodes<T>::KDTreePtInNodes(const Matrix& cloud):
		NearestNeighborSearch<T>::NearestNeighborSearch(cloud)
	{
		// build point vector and compute bounds
		BuildPoints buildPoints;
		buildPoints.reserve(cloud.cols());
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
	
	// points in nodes, priority queue
	
	template<typename T>
	KDTreePtInNodesPQ<T>::KDTreePtInNodesPQ(const Matrix& cloud):
		KDTreePtInNodes<T>::KDTreePtInNodes(cloud)
	{
	}

	template<typename T>
	typename KDTreePtInNodesPQ<T>::IndexVector KDTreePtInNodesPQ<T>::knn(const Vector& query, const Index k, const unsigned optionFlags)
	{
		typedef priority_queue<SearchElement> Queue;
		
		const bool allowSelfMatch(optionFlags & NearestNeighborSearch<T>::ALLOW_SELF_MATCH);
		
		Queue queue;
		queue.push(SearchElement(0, 0));
		IndexHeap<Index, T> heap(k);
		statistics.lastQueryVisitCount = 0;
		
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
				
				// TODO: optimise dist while going down
				const T dist(dist2<T>(node.pos, query));
				if ((dist < heap.head().value) &&
					(allowSelfMatch || (dist > numeric_limits<T>::epsilon())))
					heap.replaceHead(n, dist);
				
				// if we are at leaf, stop
				if (node.dim < 0)
					break;
				
				const T offset(query.coeff(node.dim) - node.pos.coeff(node.dim));
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
				++statistics.lastQueryVisitCount;
			}
		}
		statistics.totalVisitCount += statistics.lastQueryVisitCount;
		
		if (optionFlags & NearestNeighborSearch<T>::SORT_RESULTS)
			heap.sort();
		
		return cloudIndexesFromNodesIndexes(heap.getIndexes());
	}
	
	template struct KDTreePtInNodesPQ<float>;
	template struct KDTreePtInNodesPQ<double>;
	
	// points in nodes, stack
	
	template<typename T>
	KDTreePtInNodesStack<T>::KDTreePtInNodesStack(const Matrix& cloud):
		KDTreePtInNodes<T>::KDTreePtInNodes(cloud)
	{
	}
	
	template<typename T>
	typename KDTreePtInNodesStack<T>::IndexVector KDTreePtInNodesStack<T>::knn(const Vector& query, const Index k, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighborSearch<T>::ALLOW_SELF_MATCH);
		
		assert(nodes.size() > 0);
		assert(nodes[0].pos.size() == query.size());
		Heap heap(k);
		Vector off(Vector::Zero(nodes[0].pos.size()));
		
		statistics.lastQueryVisitCount = 0;
		
		recurseKnn(query, 0, 0, heap, off, allowSelfMatch);
		
		if (optionFlags & NearestNeighborSearch<T>::SORT_RESULTS)
			heap.sort();
		
		statistics.totalVisitCount += statistics.lastQueryVisitCount;
		// FIXME: statistics is not theard-safe
		
		return cloudIndexesFromNodesIndexes(heap.getIndexes());
	}
	
	template<typename T>
	void KDTreePtInNodesStack<T>::recurseKnn(const Vector& query, const size_t n, T rd, Heap& heap, Vector& off, const bool allowSelfMatch)
	{
		const Node& node(nodes[n]);
		const int cd(node.dim);
		
		++statistics.lastQueryVisitCount;
		
		if (cd == -2)
			return;
		
		const T dist(dist2<T>(node.pos, query));
		if ((dist < heap.head().value) &&
			(allowSelfMatch || (dist > numeric_limits<T>::epsilon()))
		)
			heap.replaceHead(n, dist);
		
		if (cd != -1)
		{
			const T old_off(off.coeff(cd));
			const T new_off(query.coeff(cd) - node.pos.coeff(cd));
			if (new_off > 0)
			{
				recurseKnn(query, childRight(n), rd, heap, off, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd < heap.head().value)
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childLeft(n), rd, heap, off, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
			else
			{
				recurseKnn(query, childLeft(n), rd, heap, off, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd < heap.head().value)
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childRight(n), rd, heap, off, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
		}
	}
	
	template struct KDTreePtInNodesStack<float>;
	template struct KDTreePtInNodesStack<double>;
	
	
	
	
	// NEW:
	
	template<typename T>
	size_t KDTreeItInLeavesStack<T>::getTreeSize(size_t elCount) const
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
	void KDTreeItInLeavesStack<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos)
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
		const size_t cutDim = argMax<T>(var);
		
		// sort
		sort(first, last, CompareDim(cutDim));
		
		// set node
		const size_t rightCount(count/2);
		const size_t leftCount(count - rightCount);
		assert(last - rightCount == first + leftCount);
		
		nodes[pos] = Node(cutDim, (first+leftCount)->pos.coeff(cutDim));
		
		//cerr << pos << " cutting on " << cutDim << " at " << (first+leftCount)->pos[cutDim] << endl;
		
		// recurse
		buildNodes(first, first + leftCount, childLeft(pos));
		buildNodes(first + leftCount, last, childRight(pos));
	}

	template<typename T>
	KDTreeItInLeavesStack<T>::KDTreeItInLeavesStack(const Matrix& cloud):
		NearestNeighborSearch<T>::NearestNeighborSearch(cloud)
	{
		// build point vector and compute bounds
		BuildPoints buildPoints;
		buildPoints.reserve(cloud.cols());
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
		//for (size_t i = 0; i < nodes.size(); ++i)
		//	cout << i << ": " << nodes[i].dim << " " << nodes[i].cutVal << endl;
	}
	
	template<typename T>
	typename KDTreeItInLeavesStack<T>::IndexVector KDTreeItInLeavesStack<T>::knn(const Vector& query, const Index k, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighborSearch<T>::ALLOW_SELF_MATCH);
		
		assert(nodes.size() > 0);
		Heap heap(k);
		Vector off(Vector::Zero(query.size()));
		
		statistics.lastQueryVisitCount = 0;
		
		recurseKnn(query, 0, 0, heap, off, allowSelfMatch);
		
		if (optionFlags & NearestNeighborSearch<T>::SORT_RESULTS)
			heap.sort();
		
		statistics.totalVisitCount += statistics.lastQueryVisitCount;
		// FIXME: statistics is not theard-safe
		
		return heap.getIndexes();
	}
	
	template<typename T>
	void KDTreeItInLeavesStack<T>::recurseKnn(const Vector& query, const size_t n, T rd, Heap& heap, Vector& off, const bool allowSelfMatch)
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
			if ((dist < heap.head().value) &&
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
				recurseKnn(query, childRight(n), rd, heap, off, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd < heap.head().value)
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childLeft(n), rd, heap, off, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
			else
			{
				recurseKnn(query, childLeft(n), rd, heap, off, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd < heap.head().value)
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childRight(n), rd, heap, off, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
		}
	}
	
	template struct KDTreeItInLeavesStack<float>;
	template struct KDTreeItInLeavesStack<double>;
}
