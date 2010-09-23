#ifndef __INDEX_HEAP_H
#define __INDEX_HEAP_H

#include "Eigen/Core"
#include "Eigen/Array"
#include <vector>
#include <algorithm>
#include <limits>

namespace Nabo
{
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
		
		IndexHeap(const Index size):
			data(size, Entry(0, std::numeric_limits<VT>::infinity()))
		{
			std::make_heap(data.begin(), data.end());
		}
		
		const Entry& head() const { return data.front(); }
		
		void replaceHead(const Index index, const Value value)
		{
			std::pop_heap(data.begin(), data.end());
			*(data.end() - 1) = Entry(index, value);
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
}

#endif // __INDEX_HEAP_H
