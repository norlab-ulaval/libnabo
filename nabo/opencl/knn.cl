T heapHeadValue(const HeapEntry* heap)
{
	return heap->value;
}

void heapHeadReplace(HeapEntry* heap, const int index, const T value, const uint K)
{
	uint i = 0;
	for (; i < K - 1; ++i)
	{
		if (heap[i + 1].value > value)
			heap[i] = heap[i + 1];
		else
			break;
	}
	heap[i].value = value;
	heap[i].index = index;
}

void heapInit(HeapEntry* heap, const uint K)
{
	for (uint i = 0; i < K; ++i)
		heap[i].value = HUGE_VALF;
}

void heapCopy(global int* dest, const HeapEntry* heap, const uint K)
{
	for (uint i = 0; i < K; ++i)
		*dest++ = heap[K-i-1].index;
}

size_t childLeft(const size_t pos) { return 2*pos + 1; }
size_t childRight(const size_t pos) { return 2*pos + 2; }

kernel void knnBruteForce(	constant Node* nodes,
						constant T* cloud,
						constant T* query,
						global int* result,
						uint K,
						T maxError,
						uint optionFlags,
						uint indexStride)
{
	HeapEntry heap[MAX_K];
	heapInit(heap, K);
	
	const size_t queryId = get_global_id(0);
	const bool allowSelfMatch = optionFlags & ALLOW_SELF_MATCH;
	constant T* q = &query[queryId * POINT_STRIDE];
	
	for (uint index = 0; index < CLOUD_POINT_COUNT; ++index)
	{
		constant T* p = &cloud[index * POINT_STRIDE];
		T dist = 0;
		for (uint i = 0; i < DIM_COUNT; ++i)
		{
			const T diff = q[i] - p[i];
			dist += diff * diff;
		}
		if (dist < heapHeadValue(heap) &&
			(allowSelfMatch || (dist > 0))) // TODO: epsilon
			heapHeadReplace(heap, index, dist, K);
	}
	heapCopy(&result[queryId * indexStride], heap, K);
}

// for cloud and result, use DirectAccessBit and stride
// preconditions:
// 		K < MAX_K
//		stack_ptr < MAX_STACK_DEPTH
kernel void knnKDTree(	constant Node* nodes,
						constant T* cloud,
						constant T* query,
						global int* result,
						uint K,
						T maxError,
						uint optionFlags,
						uint indexStride)
{
	StackEntry stack[MAX_STACK_DEPTH];
	HeapEntry heap[MAX_K];
	T off[DIM_COUNT];
	uint stackPtr = 1;

	const size_t queryId = get_global_id(0);
	const bool allowSelfMatch = optionFlags & ALLOW_SELF_MATCH;
	//const bool doSort = optionFlags & SORT_RESULTS;
	constant T* q = &query[queryId * POINT_STRIDE];
	
	heapInit(heap, K);
	
	stack[0].op = OP_BEGIN_FUNCTION;
	stack[0].n = 0;
	stack[0].rd = 0;

	while (stackPtr != 0)
	{
		--stackPtr;
		StackEntry* s = stack + stackPtr;
		const size_t n = s->n;
		constant Node* node = nodes + n;
		const int cd = node->dim;
		switch (stack[stackPtr].op)
		{
			case OP_BEGIN_FUNCTION:
			if (cd < 0)
			{
				if (cd != -1)
				{
					const int index = -(cd + 2);
					constant T* p = &cloud[index * POINT_STRIDE];
					T dist = 0;
					for (uint i = 0; i < DIM_COUNT; ++i)
					{
						const T diff = q[i] - p[i];
						dist += diff * diff;
					}
					if (dist < heapHeadValue(heap) &&
						(allowSelfMatch || (dist > 0))) // TODO: epsilon
						heapHeadReplace(heap, index, dist, K);
				}
			}
			else
			{
				s->old_off = off[cd];
				s->new_off = q[cd] - node->cutVal;
				(s+1)->op = OP_BEGIN_FUNCTION;
				s->op = OP_REC1;
				if (s->new_off > 0)
				{
					(s+1)->n = childRight(n);
					s->other_n = childLeft(n);
				}
				else
				{
					(s+1)->n = childLeft(n);
					s->other_n = childRight(n);
				}
				(s+1)->rd = s->rd;
				stackPtr+=2;
			}
			break;
			
			case OP_REC1:
			s->rd += - (s->old_off*s->old_off) + (s->new_off*s->new_off);
			if (s->rd * maxError < heapHeadValue(heap))
			{
				off[cd] = s->new_off;
				(s+1)->op = OP_BEGIN_FUNCTION;
				(s+1)->n = s->other_n;
				(s+1)->rd = s->rd;
				s->op = OP_REC2;
				stackPtr+=2;
			}
			break;
			
			case OP_REC2:
			off[cd] = s->old_off;
			break;
		}
	}
	
	heapCopy(&result[queryId * indexStride], heap, K);
}

