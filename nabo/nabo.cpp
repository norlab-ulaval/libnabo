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
	
	template struct NearestNeighborSearch<float>;
	template struct NearestNeighborSearch<double>;
}
