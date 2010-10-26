kernel void knnBruteForce(const global T* cloud,
						const global T* query,
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
	const global T* q = &query[queryId * POINT_STRIDE];
	
	for (uint index = 0; index < CLOUD_POINT_COUNT; ++index)
	{
		global T* p = &cloud[index * POINT_STRIDE];
		T dist = 0;
		for (uint i = 0; i < DIM_COUNT; ++i)
		{
			const T diff = q[i] - p[i];
			dist += diff * diff;
		}
		if (dist < heapHeadValue(heap) &&
			(allowSelfMatch || (dist > EPSILON)))
			heapHeadReplace(heap, index, dist, K);
	}
	heapCopy(&result[queryId * indexStride], heap, K);
}
