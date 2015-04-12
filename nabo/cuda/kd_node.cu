//CUDA runtime for KD_nodes



/*Roadmap: 

			(We are here)
 			 |
			\ /
			 v
Core functionality -> Dynamic Parralelism Experimentation -> Linking to existing libnabo framework -> optimization -> finalization

/Optimization |= (Search key sorting, linearizing the KD tree, improving GPU caching for node heaps)/

EST: Probably a month? Unsure. I need to look through the rest of the SDK.*/
#define maximum_depth 22
#define dim_count 3
#define K_size 16
#ifndef FLOAT_MAX 
#define FLOAT_MAX 33554430.0f
#endif 
#define BLOCK_SIZE 32
#define max_rad 256
//If coordinates are within 5% of eachother when compared to their cluster maximimum, treat them as a single point. To be used later
#define max_error 0.05f

#define OFFSIDE 0
#define ONSIDE 1
#define POINT_STRIDE 3
struct point
{
	float data[dim_count];
};
struct heap_entry
{
	float value;
	unsigned char index;
};
struct stack_entry{
	size_t n;
	uint state;
};

__device__ float heapHeadValue(heap_entry* h)
{
	return h->value;
}

__device__ heap_entry* heapHeadReplace(heap_entry* h, const int index, const float value, const uint K)
{
	uint i = 0;
	for (; i < K - 1; ++i)
	{
		if (h[i + 1].value > value)
			h[i] = h[i + 1];
		else
			break;
	}
	h[i].value = value;
	h[i].index = index;
	return h;
}

__device__ heap_entry *heapInit(const uint K)
{
	heap_entry *h;
	for (uint i = 0; i < K; ++i)
		h[i].value = FLOAT_MAX;
	return h;
}
struct kd_node
{
	//Which dimension
	unsigned int dim;
	//At what value was this node split?
	int cutVal;
	//The index of the current node
	int index;
};


#define inx_size 12
struct /*__align__(inx_size)*/ indx
{
	//The points of the KD tree
	point *pts;
	//The linked nodes
	const kd_node *nodes;
};

//Just a utility function for converting an int equal to zero to one, and vice versa. Its not well optimized, but it was quick to write :P Would be better with bitshifts
__device__ int flip(int in)
{
	return abs(in - 1);
}
__device__ unsigned int childLeft(const unsigned int pos) { return 2*pos + 1; }
__device__ unsigned int childRight(const unsigned int pos) { return 2*pos + 2; }
struct maxAB
{
	float A,B;
	int indx_a, indx_b;
};

//Clamp the value to 1 or 0
__device__ static float intensity(float a)
{
	return fmax(1,fmin(0,fabs(a)));
}
struct heap
{
	heap_entry *entries;
	int current_count;
};
//If dynamic parrallelism is not available, default to compute model 3_2. Eg: The early 700 series
#ifndef CM3_5
#define CM3_2
#endif
//Used to see if we're within bounds and are ready to jump a node
__device__ unsigned int withinBounds(int cd, point q, point p, float heapHeadVal, float maxRad, float maxError)
{
	float diff = q.data[cd] -p.data[cd];
	float side2 = diff*diff;
	if ((side2 <= maxRad) &&
		(side2 * maxError < heapHeadVal))
	{ 
		return 1;
	}
	return 0;
}
//Used for warp devices. One if distance is greater than zero. Returns 0 or 1 
__device__ unsigned int nodeMinor(int cd, point q, point p)
{
	float diff = q.data[cd] -p.data[cd];
	return (unsigned int)intensity(diff); 
	
}
//Implementation details: http://on-demand.gputechconf.com/gtc/2012/presentations/S0079-Warped-Parallel-Nearest-Neighbor-Searches-Using-KD-Trees.pdf
__device__ void recursive_warp_search(const indx static_data, const point query_point,  unsigned int _Mask, heap *output, 
					uint *stackpointer, stack_entry *stack, stack_entry *s)
{
	//Go up one
	--stackpointer;
	const size_t n = s->n;
	const kd_node node = static_data.nodes[n];
	const int cd = node.cutVal;
	//Continue doesn't do anything anymore since we're in a __device__ function (Not __global__), and there is no while loop
	/*if (cd == -2)
		continue;*/
	const int index = node.index;
	point p = static_data.pts[index];
	// compute new distance and update if lower
	float dist = 0;
	for (uint i = 0; i < dim_count; ++i)
	{
		const float diff = query_point.data[i] - p.data[i];
		dist += diff * diff;
	}
	if ((dist <= max_rad) &&
		(dist < heapHeadValue(output->entries)) &&
		(dist > (float)max_error)){
		output->entries = heapHeadReplace(output->entries, index, dist, K_size);output->current_count++;}
		// look for recursion
	//Let the warp group decide which way we want to travel next
	_Mask = _Mask & __ballot(withinBounds(cd, query_point,p,heapHeadValue(output->entries), max_rad, max_error));
	
	//This is going to be a very large amount of nested if statements o.O Luckily nvcc cleans it up a lot :D I'll eventually make it look nicer here as well

	//Did we go through the first branch?
	bool check = false;
	//If the group decided to travel off side first, continue
	if(_Mask!=0)
	{
		check = true;
		stackpointer++;
		_Mask = _Mask & __ballot(nodeMinor(cd,query_point,p));
		//Make a copy of S incase we want to branch right later on
		stack_entry s_copy1 = *s;
		s->state = OFFSIDE;
		//Branch left first
		if(_Mask!=0)
		{
			s->n = childLeft(n);	
			//Make another copy of s to come back to late
			stack_entry s_copy2 = *s;
			_Mask = _Mask & __ballot(nodeMinor(cd,query_point,p)); 
			*s = stack[*stackpointer];
			s->state = ONSIDE;	
			//We're going to branch to the right	
			stackpointer++;		
			if(_Mask !=0)
			{
				s->n = childRight(n);
				//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
				recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
				*s = s_copy2;
				//If any threads want to go to the left			
				if(__any(flip(nodeMinor(cd,query_point,p))))
				{
					//If I did not want this decision
					if(!flip(nodeMinor(cd,query_point,p)))
					{ goto exit; }
					s->n = childLeft(n);
					//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
					recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack,s);
				}
			}
			//If the group would rather branch left second
			else
			{
				s->n = childLeft(n);
				//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
				recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
				*s = s_copy2;
				//If any threads want to go to the left			
				if(__any(flip(nodeMinor(cd,query_point,p))))
				{
					//If I did not want this decision
					if(!flip(nodeMinor(cd,query_point,p)))
					{ goto exit; }
					s->n = childRight(n);
					//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
					recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack,s);
				}				
			}	
			stackpointer++;
			//We just came out of that branch, so lets set where we currently are in the stack back to s
			*s = s_copy1;
			//If any threads wanted to go right instead of left, branch that way now. In a worst case scenario, half the threads from the group are now gone
			if(__any(flip(nodeMinor(cd,query_point,p))))
			{
				//If I did not want this decision
				if(!flip(nodeMinor(cd,query_point,p)))
				{ goto exit; }
				s->n = childRight(n);
			}
		}
		//Branch right first
		else if(_Mask== 0)
		{
			s->n = childRight(n);	
			//Make another copy of s to come back to late
			stack_entry s_copy2 = *s;
			_Mask = _Mask & __ballot(nodeMinor(cd,query_point,p)); 
			*s = stack[*stackpointer];
			s->state = ONSIDE;	
			//We're going to branch to the right	
			stackpointer++;		
			if(_Mask !=0)
			{
				s->n = childRight(n);
				//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
				recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
				*s = s_copy2;
				//If any threads want to go to the left			
				if(__any(flip(nodeMinor(cd,query_point,p))))
				{
					//If I did not want this decision
					if(!flip(nodeMinor(cd,query_point,p)))
					{ goto exit; }
					s->n = childLeft(n);
					//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
					recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
				}
			}
			//If the group would rather branch left second
			else
			{
				s->n = childLeft(n);
				//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
				recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
				*s = s_copy2;
				//If any threads want to go to the left			
				if(__any(flip(nodeMinor(cd,query_point,p))))
				{
					//If I did not want this decision
					if(!flip(nodeMinor(cd,query_point,p)))
					{ goto exit; }
					s->n = childRight(n);
					//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
					recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
				}				
			}	
			//We just came out of that branch, so lets set where we currently are in the stack back to s
			*s = s_copy1;
			//If any threads wanted to go right instead of left, branch that way now. In a worst case scenario, half the threads from the group are now gone
			if(__any(flip(nodeMinor(cd,query_point,p))))
			{
				s->n = childLeft(n);
			}
		}		
	}
	//We want to branch onside and not offside
	else if((_Mask== 0) || ((check == true) && (__any(flip(nodeMinor(cd,query_point,p))))))
	{
		//Make another copy of s to come back to late
		stack_entry s_copy2 = *s;
		_Mask = _Mask & __ballot(nodeMinor(cd,query_point,p)); 
		*s = stack[*stackpointer];
		s->state = ONSIDE;	
		//We're going to branch to the right	
		stackpointer++;		
		if(_Mask !=0)
		{
			s->n = childRight(n);
			//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
			recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
			*s = s_copy2;
			//If any threads want to go to the left			
			if(__any(flip(nodeMinor(cd,query_point,p))))
			{
				//If I did not want this decision
				if(!flip(nodeMinor(cd,query_point,p)))
				{ goto exit; }
				s->n = childLeft(n);
				//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
				recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
			}
		}
		//If the group would rather branch left second
		else
		{
			s->n = childLeft(n);
			//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
			recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
			*s = s_copy2;
			//If any threads want to go to the left			
			if(__any(flip(nodeMinor(cd,query_point,p))))
			{
				//If I did not want this decision
				if(!flip(nodeMinor(cd,query_point,p)))
				{ goto exit; }				
				s->n = childRight(n);
				//Execute branch down. Eventually this will be doing through dynamic parrallelism instead of direct recursion				
				recursive_warp_search(static_data, query_point, _Mask, output, stackpointer, stack, s);
			}				
		}	
	}
exit:
	return;
	//TODO: Sort
	
} 
/*Kernel is to be executed as 32x1
indx is pre malloced and copied to the GPU to avoid memory bottlenecks. Query points is copied per iteration.
Uses a warped ballot system. Preferable for clustered points that are closely together.
Make sure the thread group size is equal to the size of the cluster & is a multiple of 32*/
__global__ void clustered_search(indx static_data, const point *query_points, int *indices,  heap *ret, int query_amt)
{
	stack_entry stack[maximum_depth];
	//Global thread number
	int thread_num = blockIdx.x * BLOCK_SIZE + threadIdx.x;
	heap myHeap;
	myHeap.entries = heapInit(K_size);
	myHeap.current_count = 0;
	//Start at root node
	stack_entry* s = stack;
	uint *startpos;
	*startpos = 0;
	recursive_warp_search(static_data, query_points[thread_num], 1, &myHeap, startpos,s,stack);
	ret[thread_num] = myHeap;
} 

