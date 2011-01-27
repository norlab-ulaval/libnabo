kernel void knnBruteForce(const global T* cloud,
						const global T* query,
						global int* indices,
						global T* dists2,
						const uint K,
						const T maxError,
						const uint optionFlags,
						const uint indexStride,
						const uint dists2Stride,
						const uint pointCount)
{
	HeapEntry heap[MAX_K];
	heapInit(heap, K);
	
	const size_t queryId = get_global_id(0);
	const bool allowSelfMatch = optionFlags & ALLOW_SELF_MATCH;
	const bool doSort = optionFlags & SORT_RESULTS;
	const global T* q = &query[queryId * POINT_STRIDE];
	
	for (uint index = 0; index < pointCount; ++index)
	{
		const global T* p = &cloud[index * POINT_STRIDE];
		T dist = 0;
		for (uint i = 0; i < DIM_COUNT; ++i)
		{
			const T diff = q[i] - p[i];
			dist += diff * diff;
		}
		if (dist < heapHeadValue(heap) &&
			(allowSelfMatch || (dist > (T)EPSILON)))
			heapHeadReplace(heap, index, dist, K);
	}
	
	if (doSort)
		heapSort(heap);
	heapCopy(&indices[queryId * indexStride], &dists2[queryId * dists2Stride], heap, K);
}
