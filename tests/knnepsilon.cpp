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
#include "helpers.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

#ifdef _POSIX_TIMERS 
namespace boost
{
	/*
		High-precision timer class, using gettimeofday().
		The interface is a subset of the one boost::timer provides,
		but the implementation is much more precise
		on systems where clock() has low precision, such as glibc.
	*/
	struct timer
	{
		typedef unsigned long long Time;
		
		timer():_start_time(curTime()){ } 
		void restart() { _start_time = curTime(); }
		double elapsed() const                  // return elapsed time in seconds
		{ return  double(curTime() - _start_time) / double(1000000000); }

	private:
		Time curTime() const {
			struct timespec ts;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
			return Time(ts.tv_sec) * Time(1000000000) + Time(ts.tv_nsec);
		}
		Time _start_time;
	};
}
#else // _POSIX_TIMERS
	#include <boost/timer.hpp>
#endif // _POSIX_TIMERS

using namespace std;
using namespace Nabo;

template<typename T>
void doTestEpsilon(const char *fileName, const int K, const int method, const int searchCount)
{
	typedef Nabo::NearestNeighbourSearch<T> NNS;
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
	
	// check methods together
	const int itCount(method >= 0 ? method : d.cols() * 2);
	
	// create big query
	const Matrix q(createQuery<T>(d, itCount, method));
	IndexMatrix indexes_bf(K, q.cols());
	Matrix dists2_bf(K, q.cols());
	
	NNS* nns = NNS::createKDTreeLinearHeap(d, std::numeric_limits<typename NNS::Index>::max(), NNS::TOUCH_STATISTICS);
	for (T epsilon = 0; epsilon < 2; epsilon += 0.01)
	{
		double duration(0);
		double touchStats(0);
		for (int s = 0; s < searchCount; ++s)
		{
			boost::timer t;
			touchStats += nns->knn(q, indexes_bf, dists2_bf, K, epsilon, NNS::ALLOW_SELF_MATCH);
			duration += t.elapsed();
		}
		cout << epsilon << " " << duration/double(searchCount) << " " << touchStats/double(searchCount) << endl;
	}
	delete nns;
}


int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		cerr << "Usage " << argv[0] << " DATA K METHOD SEARCH_COUNT" << endl;
		return 1;
	}
	
	const int K(atoi(argv[2]));
	const int method(atoi(argv[3]));
	const int searchCount(atoi(argv[4]));
	
	//cout << "Float: (epsilon, average duration)\n";
	//doTestEpsilon<float>(argv[1], K, method, searchCount);
	cout << "epsilon average_duration search_count\n";
	doTestEpsilon<double>(argv[1], K, method, searchCount);
	
	return 0;
}