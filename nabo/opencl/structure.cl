#define INVALID_NODE -1

#define OP_BEGIN_FUNCTION = 0
#define OP_REC1 = 1
#define OP_REC2 = 2

// T is the scalar type

typedef struct { 
	int dim; // -1 == invalid, <= -2 = index of pt
	T cutVal; //!< for split node, split value 
} Node;

typedef struct {
	uint op;
	size_t n;
	size_t other_n;
	T rd;
	T old_off;
	T new_off;
} StackEntry;

typedef struct {
	T value;
	int index;
} HeapEntry;