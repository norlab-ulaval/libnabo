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
		throw std::runtime_error(startMsg + "must be a multi-dimensional array");
	const int nDim = PyArray_NDIM(cloudObj);
	if (nDim != 2)
		throw std::runtime_error(startMsg + "must be a two-dimensional array");
	if (PyArray_TYPE(cloudObj) != NPY_FLOAT64)
		throw std::runtime_error(startMsg + "must hold doubles");
	if (!PyArray_CHKFLAGS(cloudObj, NPY_C_CONTIGUOUS) && !PyArray_CHKFLAGS(cloudObj, NPY_F_CONTIGUOUS))
		throw std::runtime_error(startMsg + "must be a continuous array");
}

MappedEigenDoubleMatrix* eigenFromBoostPython(const object cloudIn, const char *paramName)
{
	int dimCount, pointCount;
	const PyObject *cloudObj(cloudIn.ptr());
	
	checkPythonArray(cloudObj, paramName);
	
	matrixSizeFromPythonArray(cloudObj, dimCount, pointCount);
	double* cloudData(reinterpret_cast<double*>(PyArray_DATA(cloudObj)));
	
	return new MappedEigenDoubleMatrix(cloudData, dimCount, pointCount);
}

void eigenFromBoostPython(NNSNabo::Matrix& cloudOut, const object cloudIn, const char *paramName)
{
	int dimCount, pointCount;
	const PyObject *cloudObj(cloudIn.ptr());
	
	checkPythonArray(cloudObj, paramName);
	
	matrixSizeFromPythonArray(cloudObj, dimCount, pointCount);
	cloudOut.resize(dimCount, pointCount);
	
	memcpy(cloudOut.data(), PyArray_DATA(cloudObj), pointCount*dimCount*sizeof(double));
}

class NearestNeighbourSearch
{
public:
	NearestNeighbourSearch(const SearchType searchType, const object pycloud, const Index dim = maxI, const dict params = dict())
	{
		// build cloud
		eigenFromBoostPython(cloud, pycloud, "cloud");
		
		// build params
		Nabo::Parameters _params;
		object it = params.iteritems();
		for(int i = 0; i < len(params); ++i)
		{
			const tuple item(it.attr("next")());
			const std::string key = extract<std::string>(item[0]);
			const object val(item[1]);
			const std::string valType(val.ptr()->ob_type->tp_name);
			if (valType == "int")
			{
				const int iVal = extract<int>(val);
				if (iVal >= 0)
					_params[key] = (unsigned)iVal;
				else
					_params[key] = iVal;
			}
		}
		
		// create search
		nns = NNSNabo ::create(cloud, dim, searchType, 0, _params);
	}
	
	~NearestNeighbourSearch()
	{
		delete nns;
	}
	
	tuple knn(const object query, const Index k = 1, const double epsilon = 0, const unsigned optionFlags = 0, const double maxRadius = infD)
	{
		// map query and create output matrices
		MappedEigenDoubleMatrix* mappedQuery(eigenFromBoostPython(query, "query"));
		NNSNabo::IndexMatrix indexMatrix(k, mappedQuery->cols());
		NNSNabo::Matrix dists2Matrix(k, mappedQuery->cols());
		
		// do the search
		nns->knn(*mappedQuery, indexMatrix, dists2Matrix, k, epsilon, optionFlags, maxRadius);
		
		// build resulting python types
		npy_intp retDims[2] = { mappedQuery->cols(), k };
		const int dataCount(k * mappedQuery->cols());
		PyObject* dists2 = PyArray_EMPTY(2, retDims, PyArray_DOUBLE, PyArray_ISFORTRAN(query.ptr()));
		memcpy(PyArray_DATA(dists2), dists2Matrix.data(), dataCount*sizeof(double));
		PyObject* indices = PyArray_EMPTY(2, retDims, PyArray_INT, PyArray_ISFORTRAN(query.ptr()));
		memcpy(PyArray_DATA(indices), indexMatrix.data(), dataCount*sizeof(int));
		delete mappedQuery;
		
		// return results
		return make_tuple(object(handle<>(indices)), object(handle<>(dists2)));
	}
	
protected:
	NNSNabo *nns;
	NNSNabo::Matrix cloud;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(knn_overloads, knn, 1, 4)

BOOST_PYTHON_MODULE(pynabo)
{
	import_array();
	
	enum_<SearchType>("SearchType")
		.value("BRUTE_FORCE", NNSNabo::BRUTE_FORCE)
		.value("KDTREE_LINEAR_HEAP", NNSNabo::KDTREE_LINEAR_HEAP)
		.value("KDTREE_TREE_HEAP", NNSNabo::KDTREE_TREE_HEAP)
	;
	
	enum_<SearchOptionFlags>("SearchOptionFlags")
		.value("ALLOW_SELF_MATCH", NNSNabo::ALLOW_SELF_MATCH)
		.value("SORT_RESULTS", NNSNabo::SORT_RESULTS)
	;
	
	class_<NearestNeighbourSearch>("NearestNeighbourSearch",
		init<const SearchType, const object, optional<const Index, const dict> >())
		.def("knn", &NearestNeighbourSearch::knn, knn_overloads())
	;
}