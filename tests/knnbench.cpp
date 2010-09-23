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
typedef Nabo::KDTree<double> KDTD;

inline Vector createQuery(const Matrix& d, const KDTD& kdt, const int i, const int method)
{
	if (method == -1)
	{
		Vector q = d.col(i % d.cols());
		if (i < d.cols())
			q.cwise() += 0.01;
		else
			q.cwise() -= 0.01;
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
	const int itCount(method != -1 ? method : d.cols() * 2);
	
	// compare KDTree with brute force search
	if (K >= d.cols())
	{
		cerr << "Requested more nearest neighbor than points in the data set" << endl;
		return 2;
	}
	
	cout << "Nabo" << endl;
	cout << "\tconstruction: ";
	boost::progress_timer* t = new boost::progress_timer;
	KDTD kdt(d);
	delete t;
	
	// create queries
	Matrix q(d.rows(), itCount);
	for (int i = 0; i < itCount; ++i)
	{
		q.col(i) = createQuery(d, kdt, i, method);
	}
	
	// KDTree
	srand(0);
	cout << "\texecution: ";
	{
		boost::progress_timer t;
		for (int i = 0; i < itCount; ++i)
		{
			IndexVector indexes_kdtree(kdt.knn(q.col(i), K, 0));
		}
		
	}
	cerr << "\tstats kdtree: "
		<< kdt.getStatistics().totalVisitCount << " on "
		<< itCount * d.cols() << " ("
		<< double(100 * kdt.getStatistics().totalVisitCount) /  double(itCount * d.cols()) << " %"
		<< ")" << endl;
	
	// ANN stuff
	cout << "\nANN" << endl;
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