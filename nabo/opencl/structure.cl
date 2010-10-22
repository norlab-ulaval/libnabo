#define MAX_K 20
#define MAX_STACK_DEPTH 10

#define INVALID_NODE -1

#define OP_BEGIN_FUNCTION = 0
#define OP_REC1 = 1
#define OP_REC2 = 2

typedef float Scalar;

typedef float4 Point;

typedef struct { 
	int dim; // -1 == invalid, <= -2 = index of pt
	Scalar cutVal; //!< for split node, split value 
} Node;

typedef struct {
	uint op;
	size_t n;
	size_t other_n;
	Scalar rd;
	Scalar old_off;
	Scalar new_off;
} StackEntry;

typedef struct {
	Scalar value;
	int index;
} HeapEntry;