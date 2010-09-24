#include "nabo/nabo.h"
#include "ANN.h"
#include <boost/progress.hpp> 
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
typedef Nabo::NearestNeighborSearch<float> NNS;
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
typedef Nabo::KDTreeUnbalancedPtInLeavesImplicitBoundsStack<double> KDTD5;
typedef Nabo::KDTreeUnbalancedPtInLeavesExplicitBoundsStack<double> KDTD6;


inline Vector createQuery(const Matrix& d, const KDTD1& kdt, const int i, const int method)
{
	if (method < 0)
	{
		Vector q = d.col(i % d.cols());
		float absBound = 0;
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

template<typename T>
void doBench(const char* title, const Matrix& d, const Matrix& q, const Index K, const int itCount)
{
	boost::progress_timer* t;
	
	cout << title << endl;
	
	cout << "\tconstruction: ";
	t = new boost::progress_timer;
	T nns(d);
	delete t;
	
	cout << "\texecution: ";
	t = new boost::progress_timer;
	nns.knnM(q, K, 0);
	delete t;
	cout << "\tstats kdtree: "
		<< nns.getStatistics().totalVisitCount << " on "
		<< (long long)(itCount) * (long long)(d.cols()) << " ("
		<< (100. * double(nns.getStatistics().totalVisitCount)) /  (double(itCount) * double(d.cols())) << " %"
		<< ")\n" << endl;
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cerr << "Usage " << argv[0] << " DATA K METHOD" << endl;
		return 1;
	}
	
	const Matrix d(load<double>(argv[1]));
	const Index K(atoi(argv[2]));
	const int method(atoi(argv[3]));
	const int itCount(method >= 0 ? method : d.cols() * 2);
	
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
	
	doBench<KDTD1>("Nabo, pt in nodes, priority, balance variance", d, q, K, itCount);
	doBench<KDTD2>("Nabo, pt in nodes, stack, balance variance", d, q, K, itCount);
	doBench<KDTD3>("Nabo, balanced, stack, pt in leaves only, balance variance", d, q, K, itCount);
	doBench<KDTD4>("Nabo, balanced, stack, pt in leaves only, balance cell aspect ratio", d, q, K, itCount);
	doBench<KDTD5>("Nabo, unbalanced, stack, pt in leaves only, implicit bounds, ANN_KD_SL_MIDPT", d, q, K, itCount);
	doBench<KDTD6>("Nabo, unbalanced, points in leaves, stack, explicit bounds, ANN_KD_SL_MIDPT", d, q, K, itCount);
	
	// ANN stuff
	cout << "ANN" << endl;
	cout << "\tconstruction: ";
	boost::progress_timer* t = new boost::progress_timer;
	const int ptCount(d.cols());
	const double **pa = new const double *[d.cols()];
	for (int i = 0; i < ptCount; ++i)
		pa[i] = &d.coeff(0, i);
	ANNkd_tree* ann_kdt = new ANNkd_tree(const_cast<double**>(pa), ptCount, d.rows());
	delete t;
	
	cout << "\texecution search: ";
	{
		boost::progress_timer t;
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
	}
	cout << "\texecution pri-search: ";
	{
		boost::progress_timer t;
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
	}
	
	return 0;
}