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
typedef Nabo::KDTreePtInNodesPQ<double> KDTD1;
typedef Nabo::KDTreePtInNodesStack<double> KDTD2;
typedef Nabo::KDTreeItInLeavesStack<double> KDTD3;

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
	
	boost::progress_timer* t;
	
	// KDTD1
	cout << "Nabo priority" << endl;
	cout << "\tconstruction: ";
	t = new boost::progress_timer;
	KDTD1 kdt1(d);
	delete t;
	
	// create queries
	Matrix q(d.rows(), itCount);
	for (int i = 0; i < itCount; ++i)
	{
		q.col(i) = createQuery(d, kdt1, i, method);
	}
	
	srand(0);
	cout << "\texecution: ";
	t = new boost::progress_timer;
	kdt1.knnM(q, K, 0);
	delete t;
	cout << "\tstats kdtree: "
		<< kdt1.getStatistics().totalVisitCount << " on "
		<< itCount * d.cols() << " ("
		<< (100. * double(kdt1.getStatistics().totalVisitCount)) /  (double(itCount) * double(d.cols())) << " %"
		<< ")\n" << endl;
	
	// KDTD2
	cout << "Nabo stack" << endl;
	cout << "\tconstruction: ";
	t = new boost::progress_timer;
	KDTD2 kdt2(d);
	delete t;
	
	srand(0);
	cout << "\texecution: ";
	t = new boost::progress_timer;
	kdt2.knnM(q, K, 0);
	delete t;
	cout << "\tstats kdtree: "
		<< kdt1.getStatistics().totalVisitCount << " on "
		<< itCount * d.cols() << " ("
		<< (100. * double(kdt1.getStatistics().totalVisitCount)) /  (double(itCount) * double(d.cols())) << " %"
		<< ")\n" << endl;
	
	// KDTD3
	cout << "Nabo stack pt in leaves only" << endl;
	cout << "\tconstruction: ";
	t = new boost::progress_timer;
	KDTD3 kdt3(d);
	delete t;
	
	srand(0);
	cout << "\texecution: ";
	t = new boost::progress_timer;
	kdt3.knnM(q, K, 0);
	delete t;
	cout << "\tstats kdtree: "
		<< kdt1.getStatistics().totalVisitCount << " on "
		<< itCount * d.cols() << " ("
		<< (100. * double(kdt1.getStatistics().totalVisitCount)) /  (double(itCount) * double(d.cols())) << " %"
		<< ")\n" << endl;
	
	
	// ANN stuff
	cout << "ANN" << endl;
	cout << "\tconstruction: ";
	t = new boost::progress_timer;
	const int ptCount(d.cols());
	const double **pa = new const double *[d.cols()];
	for (int i = 0; i < ptCount; ++i)
		pa[i] = &d.coeff(0, i);
	ANNkd_tree* ann_kdt = new ANNkd_tree(const_cast<double**>(pa), ptCount, d.rows());
	delete t;
	
	srand(0);
	cout << "\texecution: ";
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