#ifndef __INDEX_HEAP_H
#define __INDEX_HEAP_H

#include "Eigen/Core"
#include "Eigen/Array"
#include <vector>
#include <algorithm>
#include <limits>

namespace Nabo
{
	// balanced-tree implementation of heap
	template<typename IT, typename VT>
	struct IndexHeap
	{
		typedef IT Index;
		typedef VT Value;
		
		struct Entry
		{
			IT index;
			VT value;
			
			Entry(const IT index, const VT value): index(index), value(value) {} 
			friend bool operator<(const Entry& e0, const Entry& e1) { return e0.value < e1.value; }
		};
		typedef std::vector<Entry> Entries;
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		
		Entries data;
		const typename Entries::iterator insertIt;
		
		IndexHeap(const Index size):
			data(size, Entry(0, std::numeric_limits<VT>::infinity())),
			insertIt(data.end() - 1)
		{
			std::make_heap(data.begin(), data.end());
		}
		
		inline void reset()
		{
			std::fill(data.begin(), data.end(), Entry(0, std::numeric_limits<VT>::infinity()));
			std::make_heap(data.begin(), data.end());
		}
		
		inline const Entry& head() const { return data.front(); }
		
		inline void replaceHead(const Index index, const Value value)
		{
			std::pop_heap(data.begin(), data.end());
			insertIt->index = index;
			insertIt->value = value;
			push_heap(data.begin(), data.end());
		}
		
		void sort()
		{
			sort_heap (data.begin(), data.end());
		}
		
		IndexVector getIndexes() const
		{
			IndexVector indexes(data.size());
			for (size_t i = 0; i < data.size(); ++i)
				indexes.coeffRef(i) = data[i].index;
			return indexes;
		}
	};
	
	// brute-force implementation of heap
	template<typename IT, typename VT>
	struct IndexHeapBruteForceVector
	{
		typedef IT Index;
		typedef VT Value;
		
		struct Entry
		{
			IT index;
			VT value;
			
			Entry(const IT index, const VT value): index(index), value(value) {} 
			friend bool operator<(const Entry& e0, const Entry& e1) { return e0.value < e1.value; }
		};
		typedef std::vector<Entry> Entries;
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		
		Entries data;
		const typename Entries::iterator headIt;
		const size_t sizeMinusOne;
		
		IndexHeapBruteForceVector(const Index size):
			data(size, Entry(0, std::numeric_limits<VT>::infinity())),
			headIt(data.end() - 1),
			sizeMinusOne(data.size() - 1)
		{
		}
		
		inline void reset()
		{
			for (typename Entries::iterator it(data.begin()); it != data.end(); ++it)
				it->value = std::numeric_limits<VT>::infinity();
		}
		
		inline const Entry& head() const { return *headIt; }
		
		inline void replaceHead(const Index index, const Value value)
		{
			register size_t i;
			for (i = sizeMinusOne; i > 0; --i)
			{
				if (data[i-1].value > value)
					data[i] = data[i-1];
				else
					break;
			}
			data[i].value = value;
			data[i].index = index;
		}
		
		void sort()
		{
			std::sort(data.begin(), data.end());
		}
		
		IndexVector getIndexes() const
		{
			IndexVector indexes(data.size());
			for (size_t i = 0; i < data.size(); ++i)
				indexes.coeffRef(i) = data[i].index;
			return indexes;
		}
	};
}

#endif // __INDEX_HEAP_H
