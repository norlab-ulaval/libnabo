#include "nabo.h"
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
		throw runtime_error(string("Cannot open file ") + fileName);
	
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
			//cout << atof(token) << " ";
			token = strtok(NULL, " \t,;"); // FIXME: non reentrant, use strtok_r
		}
		firstLine = false;
	}
	
	return Matrix::Map(&data[0], dim, data.size() / dim);
}

template<typename T>
void dumpCoordinateForSVG(const typename NearestNeighborSearch<T>::Vector coord, const float zoom = 1, const float ptSize = 1, const char* style = "stroke=\"black\" fill=\"red\"")
{
	if (coord.size() == 2)
		cout
			<< "<circle cx=\"" << zoom*coord(0)
			<< "\" cy=\"" << zoom*coord(1)
			<< "\" r=\"" << ptSize
			<< "\" stroke-width=\"" << 0.2 * ptSize
			<< "\" " << style << "/>" << endl;
	else
		assert(false);
}

int main(int argc, char* argv[])
{
	typedef Nabo::NearestNeighborSearch<float>::Matrix Matrix;
	typedef Nabo::NearestNeighborSearch<float>::Vector Vector;
	typedef Nabo::NearestNeighborSearch<float>::Index Index;
	typedef Nabo::NearestNeighborSearch<float>::Indexes Indexes;
	typedef Nabo::BruteForceSearch<float> BFSF;
	typedef Nabo::KDTree<float> KDTF;
	
	if (argc != 2)
	{
		cerr << "Usage " << argv[0] << " DATA" << endl;
		return 1;
	}
	
	Matrix d(load<float>(argv[1]));
	BFSF bfs(d);
	KDTF kdt(d);
	const Index K(10);
	
	// uncomment to compare KDTree with brute force search
	assert(K < d.size());
	const int itCount(100);
	for (int i = 0; i < itCount; ++i)
	{
		//Vector q(bfs.minBound.size());
		//for (int i = 0; i < q.size(); ++i)
		//	q(i) = bfs.minBound(i) + float(rand()) * (bfs.maxBound(i) - bfs.minBound(i)) / float(RAND_MAX);
		Vector q(d.col(rand() % d.cols()));
		q.cwise() += 0.01;
		//cerr << "bfs:\n";
		Indexes indexes_bf(bfs.knn(q, K, false));
		//cerr << "kdt:\n";
		Indexes indexes_kdtree(kdt.knn(q, K, false));
		//cout << indexes_bf.size() << " " << indexes_kdtree.size() << " " << K << endl;
		assert(indexes_bf.size() == indexes_kdtree.size());
		assert(indexes_bf.size() == K);
		//cerr << "\nquery:\n" << q << "\n\n";
		for (size_t j = 0; j < K; ++j)
		{
			Vector pbf(d.col(indexes_bf[j]));
			Vector pkdtree(d.col(indexes_kdtree[j]));
			//cerr << "index " << j << ": " << indexes_bf[j] << ", " << indexes_kdtree[j] << "\n";	
			//cerr << "point " << j << ": " << "\nbf:\n" << pbf << "\nkdtree:\n" << pkdtree << "\n\n";
			assert(dist2(pbf, pkdtree) < numeric_limits<float>::epsilon());
		}
	}
	cerr << "stats kdtree: "
		<< kdt.getStatistics().totalVisitCount << " on "
		<< itCount * d.cols() << " ("
		<< double(100 * kdt.getStatistics().totalVisitCount) /  double(itCount * d.cols()) << " %"
		<< ")" << endl;
	
	/*
	// uncomment to randomly get a point and find its minimum
	cout << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
	cout << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
	cout << "<svg width=\"100%\" height=\"100%\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl;
	srand(time(0));
	for (int i = 0; i < d.cols(); ++i)
		dumpCoordinateForSVG<float>(d.col(i), 100, 1);
	Vector q(bfs.minBound.size());
	for (int i = 0; i < q.size(); ++i)
		q(i) = bfs.minBound(i) + float(rand()) * (bfs.maxBound(i) - bfs.minBound(i)) / float(RAND_MAX);
	Indexes indexes_bf(bfs.knn(q, K, false));
	for (size_t i = 0; i < K; ++i)
		dumpCoordinateForSVG<float>(d.col(indexes_bf[i]), 100, 1,  "stroke=\"black\" fill=\"green\"");
	dumpCoordinateForSVG<float>(q, 100, 1,  "stroke=\"black\" fill=\"blue\"");
	cout << "</svg>" << endl;
	*/
	
	//cout << "Average KDTree visit count: " << double(totKDTreeVisitCount) * 100. / double(itCount * d.cols()) << " %" << endl;
	
	
	return 0;
}