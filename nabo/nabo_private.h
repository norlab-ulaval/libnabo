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
