#include "nabo/nabo.h"
#ifdef HAVE_ANN
	#include "ANN.h"
#endif // HAVE_ANN
#include <boost/timer.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace Nabo;

template<typename T>
typename NearestNeighborSearch<T>::Matrix load(const char *fileName)
{
	typedef typename NearestNeighborSearch<T>::Matrix Matrix;
	
	ifstream ifs(fileName);
	if (!ifs.good())
	{
		cerr << "Cannot open file "<< fileName << endl;
		exit(1);
	}
	
	vector<T> data;
	int dim(0);
	bool firstLine(true);
	
	while (!ifs.eof())
	{
		char line[1024];
		ifs.getline(line, sizeof(line));
		line[sizeof(line)-1] = 0;
		
		char *token = strtok(line, " \t,;");
		while (token)
		{
			if (firstLine)
				++dim;
			data.push_back(atof(token));
			token = strtok(NULL, " \t,;"); // FIXME: non reentrant, use strtok_r
		}
		firstLine = false;
	}
	
	return Matrix::Map(&data[0], dim, data.size() / dim);
}

typedef Nabo::NearestNeighborSearch<double>::Matrix Matrix;
typedef Nabo::NearestNeighborSearch<double>::Vector Vector;
typedef Nabo::NearestNeighborSearch<double>::Index Index;
typedef Nabo::NearestNeighborSearch<double>::IndexVector IndexVector;
typedef Nabo::NearestNeighborSearch<double> NNS;
typedef Nabo::BruteForceSearch<double> BFSD;
typedef Nabo::KDTreeBalancedPtInNodesPQ<double> KDTD1;
typedef Nabo::KDTreeBalancedPtInNodesStack<double> KDTD2;
struct KDTD3: public Nabo::KDTreeBalancedPtInLeavesStack<double>
{
	KDTD3(const Matrix& cloud):
		Nabo::KDTreeBalancedPtInLeavesStack<double>(cloud, true)
	{}
};
struct KDTD4: public Nabo::KDTreeBalancedPtInLeavesStack<double>
{
	KDTD4(const Matrix& cloud):
		Nabo::KDTreeBalancedPtInLeavesStack<double>(cloud, false)
	{}
};
typedef Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStack<double,IndexHeapSTL<int,double>> KDTD5A;
typedef Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStack<double,IndexHeapBruteForceVector<int,double>> KDTD5B;

//typedef Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapBruteForceVector<int,double>> KDTD5OB;
//typedef Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapSTL<int,double>> KDTD5OA;

template<typename T, typename Heap, int C>
struct KDTD5O: public Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T,Heap>
{
	KDTD5O(const Matrix& cloud):
		Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T,Heap>(cloud, C)
	{}
};

typedef KDTD5O<double,IndexHeapBruteForceVector<int,double>,1> KDTD5OB1;
typedef KDTD5O<double,IndexHeapBruteForceVector<int,double>,2> KDTD5OB2;
typedef KDTD5O<double,IndexHeapBruteForceVector<int,double>,2> KDTD5OB4;
typedef KDTD5O<double,IndexHeapSTL<int,double>,1> KDTD5OS1;
typedef KDTD5O<double,IndexHeapSTL<int,double>,2> KDTD5OS2;
typedef KDTD5O<double,IndexHeapSTL<int,double>,4> KDTD5OS4;

typedef Nabo::KDTreeUnbalancedPtInLeavesExplicitBoundsStack<double> KDTD6;


inline Vector createQuery(const Matrix& d, const KDTD1& kdt, const int i, const int method)
{
	if (method < 0)
	{
		Vector q = d.col(i % d.cols());
		double absBound = 0;
		for (int j = 0; j < q.size(); ++j)
			absBound += kdt.maxBound(j) - kdt.minBound(j);
		absBound /= 3 * (-method); // divided by -method
		if (i < d.cols())
			q.cwise() += absBound;
		else
			q.cwise() -= absBound;
		return q;
	}
	else
	{
		Vector q(kdt.minBound.size());
		for (int j = 0; j < q.size(); ++j)
			q(j) = kdt.minBound(j) + double(rand()) * (kdt.maxBound(j) - kdt.minBound(j)) / double(RAND_MAX);
		return q;
	}
}

struct BenchResult
{
	double creationDuration;
	double executionDuration;
	double visitCount;
	double totalCount;
	
	BenchResult():
		creationDuration(0),
		executionDuration(0),
		visitCount(0),
		totalCount(0)
	{}
	
	void operator +=(const BenchResult& that)
	{
		creationDuration += that.creationDuration;
		executionDuration += that.executionDuration;
		visitCount += that.visitCount;
		totalCount += that.totalCount;
	}
	
	void operator /=(const double factor)
	{
		creationDuration /= factor;
		executionDuration /= factor;
		visitCount /= factor;
		totalCount /= factor;
	}
};
typedef vector<BenchResult> BenchResults;

template<typename T>
BenchResult doBench(const Matrix& d, const Matrix& q, const Index K, const int itCount)
{
	BenchResult result;
	boost::timer t;
	T nns(d);
	result.creationDuration = t.elapsed();
	
	t.restart();
	nns.knnM(q, K, 0, 0);
	result.executionDuration = t.elapsed();
	
	result.visitCount = double(nns.getStatistics().totalVisitCount);
	result.totalCount = double(itCount) * double(d.cols());
	
	return result;
}

#ifdef HAVE_ANN

BenchResult doBenchANNStack(const Matrix& d, const Matrix& q, const Index K, const int itCount)
{
	BenchResult result;
	boost::timer t;
	const int ptCount(d.cols());
	const double **pa = new const double *[d.cols()];
	for (int i = 0; i < ptCount; ++i)
		pa[i] = &d.coeff(0, i);
	ANNkd_tree* ann_kdt = new ANNkd_tree(const_cast<double**>(pa), ptCount, d.rows());
	result.creationDuration = t.elapsed();
	
	t.restart();
	ANNidx nnIdx[K];
	ANNdist dists[K];
	for (int i = 0; i < itCount; ++i)
	{
		const Vector& tq(q.col(i));
		ANNpoint queryPt(const_cast<double*>(&tq.coeff(0)));
		ann_kdt->annkSearch(		// search
						queryPt,	// query point
						K,			// number of near neighbors
						nnIdx,		// nearest neighbors (returned)
						dists,		// distance (returned)
						0);			// error bound
	}
	result.executionDuration = t.elapsed();
	delete ann_kdt;
	
	return result;
}

BenchResult doBenchANNPriority(const Matrix& d, const Matrix& q, const Index K, const int itCount)
{
	BenchResult result;
	boost::timer t;
	const int ptCount(d.cols());
	const double **pa = new const double *[d.cols()];
	for (int i = 0; i < ptCount; ++i)
		pa[i] = &d.coeff(0, i);
	ANNkd_tree* ann_kdt = new ANNkd_tree(const_cast<double**>(pa), ptCount, d.rows());
	result.creationDuration = t.elapsed();
	
	t.restart();
	ANNidx nnIdx[K];
	ANNdist dists[K];
	for (int i = 0; i < itCount; ++i)
	{
		const Vector& tq(q.col(i));
		ANNpoint queryPt(const_cast<double*>(&tq.coeff(0)));
		ann_kdt->annkPriSearch(		// search
						queryPt,	// query point
						K,			// number of near neighbors
						nnIdx,		// nearest neighbors (returned)
						dists,		// distance (returned)
						0);			// error bound
	}
	result.executionDuration = t.elapsed();
	
	return result;
}

#endif // HAVE_ANN


int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		cerr << "Usage " << argv[0] << " DATA K METHOD RUN_COUNT" << endl;
		return 1;
	}
	
	const Matrix d(load<double>(argv[1]));
	const Index K(atoi(argv[2]));
	const int method(atoi(argv[3]));
	const int itCount(method >= 0 ? method : d.cols() * 2);
	const int runCount(atoi(argv[4]));
	
	// compare KDTree with brute force search
	if (K >= d.cols())
	{
		cerr << "Requested more nearest neighbor than points in the data set" << endl;
		return 2;
	}
	
	// create queries
	Matrix q(d.rows(), itCount);
	{
		// temp kdtree to get bounds
		KDTD1 kdtt(d);
		for (int i = 0; i < itCount; ++i)
		{
			q.col(i) = createQuery(d, kdtt, i, method);
		}
	}
	
	#define SELF_BENCH_COUNT 6
	#ifdef HAVE_ANN
		#define BENCH_COUNT (SELF_BENCH_COUNT+1)
	#else // HAVE_ANN
		#define BENCH_COUNT SELF_BENCH_COUNT
	#endif // HAVE_ANN
	const size_t benchCount(BENCH_COUNT);
	const char* benchLabels[benchCount] =
	{
		//doBench<KDTD1>("Nabo, pt in nodes, priority, balance variance",
		//doBench<KDTD2>("Nabo, pt in nodes, stack, balance variance",
		//doBench<KDTD3>("Nabo, balanced, stack, pt in leaves only, balance variance",
		//"Nabo, balanced, stack, pt in leaves only, balance cell aspect ratio",
		//"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, STL heap",
		//"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, brute-force vector heap",
		"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, brute-force vector heap, opt, 1 thread",
		"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, brute-force vector heap, opt, 2 threads",
		"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, brute-force vector heap, opt, 4 threads",
		"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, STL heap, opt, 1 thread",
		"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, STL heap, opt, 2 threads",
		"Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT, STL heap, opt, 4 threads",
		//"Nabo, unbalanced, points in leaves, stack, explicit bounds, ANN_KD_SL_MIDPT",
		#ifdef HAVE_ANN
		"ANN stack",
		//"ANN priority",
		#endif // HAVE_ANN
	};
	
	// do bench themselves, accumulate over several times
	BenchResults results(benchCount);
	for (int run = 0; run < runCount; ++run)
	{
		size_t i = 0;
		//results.at(i++) += doBench<KDTD1>(d, q, K, itCount);
		//results.at(i++) += doBench<KDTD2>(d, q, K, itCount);
		//results.at(i++) += doBench<KDTD3>(d, q, K, itCount);
		//results.at(i++) += doBench<KDTD4>(d, q, K, itCount);
		//results.at(i++) += doBench<KDTD5A>(d, q, K, itCount);
		//results.at(i++) += doBench<KDTD5B>(d, q, K, itCount);
		results.at(i++) += doBench<KDTD5OB1>(d, q, K, itCount);
		results.at(i++) += doBench<KDTD5OB2>(d, q, K, itCount);
		results.at(i++) += doBench<KDTD5OB4>(d, q, K, itCount);
		results.at(i++) += doBench<KDTD5OS1>(d, q, K, itCount);
		results.at(i++) += doBench<KDTD5OS2>(d, q, K, itCount);
		results.at(i++) += doBench<KDTD5OS4>(d, q, K, itCount);
		//results.at(i++) += doBench<KDTD6>(d, q, K, itCount);
		#ifdef HAVE_ANN
		results.at(i++) += doBenchANNStack(d, q, K, itCount);
		//results.at(i++) += doBenchANNPriority(d, q, K, itCount);
		#endif // HAVE_ANN
	}
	
	// print results
	cout << "Showing average over " << runCount << " runs\n\n";
	for (size_t i = 0; i < benchCount; ++i)
	{
		results[i] /= double(runCount);
		cout << "Method " << benchLabels[i] << ":\n";
		cout << "  creation duration: " << results[i].creationDuration << "\n";
		cout << "  execution duration: " << results[i].executionDuration << "\n";
		if (results[i].totalCount != 0)
		{
			cout << "  visit count: " << results[i].visitCount << "\n";
			cout << "  total count: " << results[i].totalCount << "\n";
			cout << "  precentage visit: " << (results[i].visitCount * 100.) / results[i].totalCount << "\n";
		}
		else
			cout << "  no stats for visits\n";
		cout << endl;
	}
	
	return 0;
}