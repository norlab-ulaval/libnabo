#include "nabo.h"
#include "index_heap.h"

namespace Nabo
{
	using namespace std;
	
	template<typename T>
	BruteForceSearch<T>::BruteForceSearch(const Matrix& cloud):
		NearestNeighborSearch<T>::NearestNeighborSearch(cloud)
	{
		// compute bounds
		for (size_t i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			const_cast<Vector&>(this->minBound) = this->minBound.cwise().min(v);
			const_cast<Vector&>(this->maxBound) = this->maxBound.cwise().max(v);
		}
	}
	

	template<typename T>
	typename NearestNeighborSearch<T>::IndexVector BruteForceSearch<T>::knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighborSearch<T>::ALLOW_SELF_MATCH);
		
		IndexHeap<Index, T> heap(k);
		
		for (size_t i = 0; i < this->cloud.cols(); ++i)
		{
			const T dist(dist2<T>(this->cloud.col(i), query));
			if ((dist < heap.head().value) &&
				(allowSelfMatch || (dist > numeric_limits<T>::epsilon())))
				heap.replaceHead(i, dist);
		}
		
		this->statistics.lastQueryVisitCount = this->cloud.cols();
		this->statistics.totalVisitCount += this->statistics.lastQueryVisitCount;
		
		if (optionFlags & NearestNeighborSearch<T>::SORT_RESULTS)
			heap.sort();
		
		return heap.getIndexes();
	}
	
	template struct BruteForceSearch<float>;
	template struct BruteForceSearch<double>;
}
