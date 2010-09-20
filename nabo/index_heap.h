#ifndef __INDEX_HEAP_H
#define __INDEX_HEAP_H

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
		typedef std::vector<Index> Indexes;
		
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
		
		Indexes getIndexes() const
		{
			Indexes indexes;
			indexes.reserve(data.size());
			for (typename Entries::const_iterator it(data.begin()); it != data.end(); ++it)
				indexes.push_back(it->index);
			return indexes;
		}
	};
}

#endif // __INDEX_HEAP_H
