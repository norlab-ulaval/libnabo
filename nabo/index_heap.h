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
	struct IndexHeapSTL
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
		const VT& headValueRef;
		const typename Entries::iterator insertIt;
		
		IndexHeapSTL(const size_t size):
			data(size, Entry(0, std::numeric_limits<VT>::infinity())),
			headValueRef(data.begin()->value),
			insertIt(data.end() - 1)
		{
			std::make_heap(data.begin(), data.end());
		}
		
		inline void reset()
		{
			std::fill(data.begin(), data.end(), Entry(0, std::numeric_limits<VT>::infinity()));
			std::make_heap(data.begin(), data.end());
		}
		
		inline const VT& headValue() const { return headValueRef; }
		
		inline void replaceHead(const Index index, const Value value)
		{
			std::pop_heap(data.begin(), data.end());
			insertIt->index = index;
			insertIt->value = value;
			push_heap(data.begin(), data.end());
		}
		
		inline void sort()
		{
			sort_heap (data.begin(), data.end());
		}
		
		inline IndexVector getIndexes() const
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
		const VT& headValueRef;
		const size_t sizeMinusOne;
		
		IndexHeapBruteForceVector(const size_t size):
			data(size, Entry(0, std::numeric_limits<VT>::infinity())),
			headValueRef((data.end() - 1)->value),
			sizeMinusOne(data.size() - 1)
		{
		}
		
		inline void reset()
		{
			for (typename Entries::iterator it(data.begin()); it != data.end(); ++it)
				it->value = std::numeric_limits<VT>::infinity();
		}
		
		inline const VT& headValue() const { return headValueRef; }
		
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
		
		inline void sort()
		{
			std::sort(data.begin(), data.end());
		}
		
		inline IndexVector getIndexes() const
		{
			IndexVector indexes(data.size());
			for (size_t i = 0; i < data.size(); ++i)
				indexes.coeffRef(i) = data[i].index;
			return indexes;
		}
	};
}

#endif // __INDEX_HEAP_H
