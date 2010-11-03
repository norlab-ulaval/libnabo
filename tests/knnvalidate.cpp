/*

Copyright (c) 2010, Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the author at <stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "nabo/nabo.h"
//#include "experimental/nabo_experimental.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace Nabo;

template<typename T>
typename NearestNeighbourSearch<T>::Matrix load(const char *fileName)
{
	typedef typename NearestNeighbourSearch<T>::Matrix Matrix;
	
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

template<typename T, typename NNS>
typename Nabo::NearestNeighbourSearch<T>::Vector createQuery(const typename Nabo::NearestNeighbourSearch<T>::Matrix& d, const NNS& nns, const int i, const int method)
{
	if (method < 0)
	{
		typename Nabo::NearestNeighbourSearch<T>::Vector q = d.col(i % d.cols());
		double absBound = 0;
		for (int j = 0; j < q.size(); ++j)
			absBound += nns.maxBound(j) - nns.minBound(j);
		absBound /= 3 * (-method); // divided by -method
		if (i < d.cols())
			q.cwise() += absBound;
		else
			q.cwise() -= absBound;
		return q;
	}
	else
	{
		typename Nabo::NearestNeighbourSearch<T>::Vector q(nns.minBound.size());
		for (int j = 0; j < q.size(); ++j)
			q(j) = nns.minBound(j) + double(rand()) * (nns.maxBound(j) - nns.minBound(j)) / double(RAND_MAX);
		return q;
	}
}

template<typename T>
typename NearestNeighbourSearch<T>::Matrix createQuery(const typename NearestNeighbourSearch<T>::Matrix& d, const int itCount, const int method)
{
	typedef typename NearestNeighbourSearch<T>::Matrix MatrixT;
	typedef Nabo::NearestNeighbourSearch<T> NNS;
	NNS* nns = NNS::create(d, d.rows(), typename NNS::SearchType(0));
	MatrixT q(d.rows(), itCount);
	for (int i = 0; i < itCount; ++i)
		q.col(i) = createQuery<T>(d, *nns, i, method);
	delete nns;
	return q;
}


template<typename T>
void validate(const char *fileName, const int K, const int method)
{
	typedef Nabo::NearestNeighbourSearch<T> NNS;
	typedef vector<NNS*> NNSV;
	typedef typename NNS::Matrix Matrix;
	typedef typename NNS::Vector Vector;
	typedef typename NNS::Index Index;
	typedef typename NNS::IndexVector IndexVector;
	typedef typename NNS::IndexMatrix IndexMatrix;
	
	// check if file is ok
	const Matrix d(load<T>(fileName));
	if (K >= d.cols())
	{
		cerr << "Requested more nearest neighbour than points in the data set" << endl;
		exit(2);
	}
	
	// create different methods
	NNSV nnss;
	unsigned searchTypeCount(NNS::SEARCH_TYPE_COUNT);
	#ifndef HAVE_OPENCL
	searchTypeCount -= 2;
	#endif // HAVE_OPENCL
	for (unsigned i = 0; i < searchTypeCount; ++i)
		nnss.push_back(NNS::create(d, d.rows(), typename NNS::SearchType(i)));
	//nnss.push_back(new KDTreeBalancedPtInLeavesStack<T>(d, false));
	
	
	// check methods together
	const int itCount(method != -1 ? method : d.cols() * 2);
	
	/*
	// element-by-element search
	for (int i = 0; i < itCount; ++i)
	{
		const Vector q(createQuery<T>(d, *nnss[0], i, method));
		const IndexVector indexes_bf(nnss[0]->knn(q, K, 0, NNS::SORT_RESULTS));
		for (size_t j = 1; j < nnss.size(); ++j)
		{
			const IndexVector indexes_kdtree(nnss[j]->knn(q, K, 0, NNS::SORT_RESULTS));
			if (indexes_bf.size() != K)
			{
				cerr << "Different number of points found between brute force and request" << endl;
				exit(3);
			}
			if (indexes_bf.size() != indexes_kdtree.size())
			{
				cerr << "Different number of points found between brute force and NNS type "<< j  << endl;
				exit(3);
			}
			for (size_t k = 0; k < size_t(K); ++k)
			{
				Vector pbf(d.col(indexes_bf[k]));
				//cerr << indexes_kdtree[k] << endl;
				Vector pkdtree(d.col(indexes_kdtree[k]));
				if (fabsf((pbf-q).squaredNorm() - (pkdtree-q).squaredNorm()) >= numeric_limits<float>::epsilon())
				{
					cerr << "Method " << j << ", cloud point " << i << ", neighbour " << k << " of " << K << " is different between bf and kdtree (dist " << (pbf-pkdtree).norm() << ")\n";
					cerr << "* query:\n";
					cerr << q << "\n";
					cerr << "* indexes " << indexes_bf[k] << " (bf) vs " <<  indexes_kdtree[k] << " (kdtree)\n";
					cerr << "* coordinates:\n";
					cerr << "bf: (dist " << (pbf-q).norm() << ")\n";
					cerr << pbf << "\n";
					cerr << "kdtree (dist " << (pkdtree-q).norm() << ")\n";
					cerr << pkdtree << endl;
					exit(4);
				}
			}
		}
	}
	*/
	// create big query
	// check all-in-one query
	Matrix q(createQuery<T>(d, itCount, method));
	IndexMatrix indexes_bf(K, q.cols());
	Matrix dists2_bf(K, q.cols());
	nnss[0]->knn(q, indexes_bf, dists2_bf, K, 0, NNS::SORT_RESULTS);
	assert(indexes_bf.cols() == q.cols());
	for (size_t j = 1; j < nnss.size(); ++j)
	{
		IndexMatrix indexes_kdtree(K, q.cols());
		Matrix dists2_kdtree(K, q.cols());
		nnss[j]->knn(q, indexes_kdtree, dists2_kdtree, K, 0, NNS::SORT_RESULTS);
		if (indexes_bf.rows() != K)
		{
			cerr << "Different number of points found between brute force and request" << endl;
			exit(3);
		}
		if (indexes_bf.cols() != indexes_kdtree.cols())
		{
			cerr << "Different number of points found between brute force and NNS type "<< j  << endl;
			exit(3);
		}
		
		for (int i = 0; i < q.cols(); ++i)
		{
			for (size_t k = 0; k < size_t(K); ++k)
			{
				const int pbfi(indexes_bf(k,i));
				const Vector pbf(d.col(pbfi));
				const int pkdt(indexes_kdtree(k,i));
				if (pkdt < 0 || pkdt >= d.cols())
				{
					cerr << "Method " << j << ", query point " << i << ", neighbour " << k << " of " << K << " has invalid index " << pkdt << " out of range [0:" << d.cols() << "[" << endl;
					exit(4);
				}
				const Vector pkdtree(d.col(pkdt));
				const Vector pq(q.col(i));
				const float distDiff(fabsf((pbf-pq).squaredNorm() - (pkdtree-pq).squaredNorm()));
				if (distDiff > numeric_limits<float>::epsilon())
				{
					cerr << "Method " << j << ", query point " << i << ", neighbour " << k << " of " << K << " is different between bf and kdtree (dist2 " << distDiff << ")\n";
					cerr << "* query point:\n";
					cerr << pq << "\n";
					cerr << "* indexes " << pbfi << " (bf) vs " << pkdt << " (kdtree)\n";
					cerr << "* coordinates:\n";
					cerr << "bf: (dist " << (pbf-pq).norm() << ")\n";
					cerr << pbf << "\n";
					cerr << "kdtree (dist " << (pkdtree-pq).norm() << ")\n";
					cerr << pkdtree << endl;
					cerr << "* bf neighbours:\n";
					for (int l = 0; l < K; ++l)
						cerr << indexes_bf(l,i) << " (dist " << (d.col(indexes_bf(l,i))-pq).norm() << ")\n";
					cerr << "* kdtree neighbours:\n";
					for (int l = 0; l < K; ++l)
						cerr << indexes_kdtree(l,i) << " (dist " << (d.col(indexes_kdtree(l,i))-pq).norm() << ")\n";
					exit(4);
				}
			}
		}
	}
	
// 	cout << "\tstats kdtree: "
// 		<< kdt.getStatistics().totalVisitCount << " on "
// 		<< (long long)(itCount) * (long long)(d.cols()) << " ("
// 		<< (100. * double(kdt.getStatistics().totalVisitCount)) /  (double(itCount) * double(d.cols())) << " %"
// 		<< ")\n" << endl;
	
	// delete searches
	for (typename NNSV::iterator it(nnss.begin()); it != nnss.end(); ++it)
		delete (*it);
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cerr << "Usage " << argv[0] << " DATA K METHOD" << endl;
		return 1;
	}
	
	const int K(atoi(argv[2]));
	const int method(atoi(argv[3]));
	
	validate<float>(argv[1], K, method);
	//validate<double>(argv[1], K, method);
	
	return 0;
}