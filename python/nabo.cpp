#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include "../nabo/nabo.h"
#include <iostream>
#include <cassert>

using namespace boost::python;

typedef Nabo::NNSearchD NNSNabo;
typedef NNSNabo::Index Index;
typedef NNSNabo::SearchType SearchType;
typedef NNSNabo::SearchOptionFlags SearchOptionFlags;
typedef Eigen::Map<NNSNabo::Matrix> MappedEigenDoubleMatrix;
typedef Eigen::Map<NNSNabo::IndexMatrix> MappedEigenIndexMatrix;

static const double infD = std::numeric_limits<double>::infinity();
static const Index maxI = std::numeric_limits<Index>::max();

#if PY_MAJOR_VERSION >= 3
int
init_numpy()
{
	import_array();
	return 0;
}
#else
void
init_numpy()
{
	import_array();
}
#endif

void matrixSizeFromPythonArray(const PyObject* cloudObj, int& rowCount, int& colCount)
{
	assert(PyArray_CHKFLAGS(cloudObj, NPY_C_CONTIGUOUS) || PyArray_CHKFLAGS(cloudObj, NPY_F_CONTIGUOUS));
	assert(PyArray_NDIM(cloudObj) == 2);
	const npy_intp *shape = PyArray_DIMS(cloudObj);
	if (PyArray_CHKFLAGS(cloudObj, NPY_F_CONTIGUOUS))
	{
		colCount = shape[1];
		rowCount = shape[0];
	}
	else
	{
		colCount = shape[0];
		rowCount = shape[1];
	}
}

void checkPythonArray(const PyObject* cloudObj, const char *paramName)
{
	std::string startMsg("Argument \"");
	startMsg += paramName;
	startMsg += "\" ";
	
	if (!PyArray_Check(cloudObj))
		throw s{
	init_numpy();
	
	enum_<SearchType>("SearchType", "Type of algorithm used for search.")
		.value("BRUTE_FORCE", NNSNabo::BRUTE_FORCE)
		.value("KDTREE_LINEAR_HEAP", NNSNabo::KDTREE_LINEAR_HEAP)
		.value("KDTREE_TREE_HEAP", NNSNabo::KDTREE_TREE_HEAP)
	;
	
	enum_<SearchOptionFlags>("SearchOptionFlags", "Flags you can OR when creating search.")
		.value("ALLOW_SELF_MATCH", NNSNabo::ALLOW_SELF_MATCH)
		.value("SORT_RESULTS", NNSNabo::SORT_RESULTS)
	;
	
	class_<NearestNeighbourSearch>(
		"NearestNeighbourSearch",
		"Nearest-neighbour search object, containing the data, on which you can do the knn(...) query.\n\n"
		"The data and query must be continuous numpy arrays.\n"
		"As numpy proposes both C and Fortran data orders, pynabo\n"
		"will always consider the contiguous dimension to be coordinates\n"
		"of points, regardless of order, as this provides the fastest\n"
		"possible execution. The return values of knn(...) will have\n"
		"the same order as the query and will have the different results\n"
		"of each point in the contiguous dimension."
		,
		init<const object, optional<const SearchType, const Index, const dict> >(
			"Create a nearest-neighbour search.\n\n"
			"Arguments:\n"
			"    data -- data-point cloud in which to search, must be a numpy array\n"
			"    searchType -- type of search, default: KDTREE_LINEAR_HEAP\n"
			"    dim -- number of dimensions to consider, must be lower or equal to cloud.rows(), default: dim of data\n"
			"    creationOptionFlags -- creation options, a bitwise OR of elements of CreationOptionFlags",
			args("self", "data", "searchType", "dim", "creationOptionFlags, default: 0")
		)
	)
	.def("knn", &NearestNeighbourSearch::knn,
		knn_overloads(
			args("self", "query", "k", "epsilon", "optionFlags", "maxRadius"),
			"Find the k nearest neighbours of query in data.\n\n"
			"Arguments:\n"
			"    query -- query points, must be a numpy array\n"
			"    k -- number of nearest neighbour requested, default: 1\n"
			"    epsilon -- maximal ratio of error for approximate search, 0 for exact search; has no effect if the number of neighbour found is smaller than the number requested; default: 0.\n"
			"    optionFlags -- search options, a bitwise OR of elements of SearchOptionFlags, default: 0\n"
			"    maxRadius -- maximum radius in which to search, can be used to prune search, is not affected by epsilon, default: inf\n\n"
			"Returns:\n"
			"    A tuple of two 2D numpy arrays, the first containing indices to points in data, the other containing squared distances."
		)
	)
	;
}
