#include "nabo/nabo.h"
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

int main(int argc, char* argv[])
{
	typedef Nabo::NearestNeighborSearch<float>::Matrix Matrix;
	typedef Nabo::NearestNeighborSearch<float>::Vector Vector;
	typedef Nabo::NearestNeighborSearch<float>::Index Index;
	typedef Nabo::NearestNeighborSearch<float>::Indexes Indexes;
	typedef Nabo::BruteForceSearch<float> BFSF;
	typedef Nabo::KDTree<float> KDTF;
	
	if (argc != 4)
	{
		cerr << "Usage " << argv[0] << " DATA K METHOD" << endl;
		return 1;
	}
	
	const Matrix d(load<float>(argv[1]));
	const Index K(atoi(argv[2]));
	const int method(atoi(argv[3]));
	BFSF bfs(d);
	KDTF kdt(d);
	
	// compare KDTree with brute force search
	if (K >= d.cols())
	{
		cerr << "Requested more nearest neighbor than points in the data set" << endl;
		return 2;
	}
	
	const int itCount(method != -1 ? method : d.cols() * 2);
	for (int i = 0; i < itCount; ++i)
	{
		Vector q(bfs.minBound.size());
		if (method == -1)
		{
			q = d.col(i % d.cols());
			if (i < itCount / 2)
				q.cwise() += 0.01;
			else
				q.cwise() -= 0.01;
		}
		else
		{
			for (int j = 0; j < q.size(); ++j)
				q(j) = bfs.minBound(j) + float(rand()) * (bfs.maxBound(j) - bfs.minBound(j)) / float(RAND_MAX);
		}
		Indexes indexes_bf(bfs.knn(q, K, false));
		Indexes indexes_kdtree(kdt.knn(q, K, false));
		if (indexes_bf.size() != indexes_kdtree.size())
		{
			cerr << "Different number of points found between search methods" << endl;
			return 3;
		}
		if (indexes_bf.size() != K)
		{
			cerr << "Different number of points found between brute force and request" << endl;
			return 3;
		}
		for (size_t j = 0; j < K; ++j)
		{
			Vector pbf(d.col(indexes_bf[j]));
			Vector pkdtree(d.col(indexes_kdtree[j]));
			if (fabsf(dist2(pbf, q) - dist2(pkdtree, q)) >= numeric_limits<float>::epsilon())
			{
				cerr << "Point " << j << " of " << K << " is different between bf and kdtree (dist " << dist2(pbf, pkdtree) << ")\n";
				cerr << "* query:\n";
				cerr << q << "\n";
				cerr << "* indexes " << indexes_bf[j] << " (bf) vs " <<  indexes_kdtree[j] << " (kdtree)\n";
				cerr << "* coordinates:\n";
				cerr << "bf: (dist " << dist2(pbf, q) << ")\n";
				cerr << pbf << "\n";
				cerr << "kdtree (dist " << dist2(pkdtree, q) << ")\n";
				cerr << pkdtree << endl;
				return 4;
			}
		}
	}
	cout << "stats kdtree: "
		<< kdt.getStatistics().totalVisitCount << " on "
		<< itCount * d.cols() << " ("
		<< double(100 * kdt.getStatistics().totalVisitCount) /  double(itCount * d.cols()) << " %"
		<< ")" << endl;
	
	return 0;
}