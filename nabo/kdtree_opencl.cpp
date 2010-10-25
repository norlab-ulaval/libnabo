/*
 * 
 C opyr*ight (c) 2010, Stephane Magnenat, ASL, ETHZ, Switzerland
 You can contact the author at <stephane at magnenat dot net>
 
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * Neither the name of the <organization> nor the
 * names of its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * 
 T HIS *SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
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

#ifdef HAVE_OPENCL

#include "nabo_private.h"
#include "index_heap.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <queue>
#include <algorithm>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>

/*!	\file kdtree_opencl.cpp
 \ *brief kd-tree search, opencl implementation
 \ingroup private
 */

namespace cl
{
	typedef std::vector<Device> Devices;
}

namespace Nabo
{
	//! \ingroup private
	//@{
	
	//! Maximum number of points acceptable in a query
	#define MAX_K 32
	
	using namespace std;
	
	template<typename T> struct TypeName {};
	#define DEF_TYPE_MAP(CT, CLT) \
		template<> struct TypeName<CT> { static const char name[]; }; \
		const char TypeName<CT>::name[] = CLT;
	DEF_TYPE_MAP(float, "float");
	DEF_TYPE_MAP(double, "double");
	
	
	template<typename T>
	size_t argMax(const typename NearestNeighbourSearch<T>::Vector& v)
	{
		T maxVal(0);
		size_t maxIdx(0);
		for (int i = 0; i < v.size(); ++i)
		{
			if (v[i] > maxVal)
			{
				maxVal = v[i];
				maxIdx = i;
			}
		}
		return maxIdx;
	}

	template<typename T>
	size_t KDTreeBalancedPtInLeavesStackOpenCL<T>::getTreeSize(size_t elCount) const
	{
		// FIXME: 64 bits safe stuff, only work for 2^32 elements right now
		assert(elCount > 0);
		elCount --;
		size_t count = 0;
		int i = 31;
		for (; i >= 0; --i)
		{
			if (elCount & (1 << i))
				break;
		}
		for (int j = 0; j <= i; ++j)
			count |= (1 << j);
		count <<= 1;
		count |= 1;
		return count;
	}
	
	template<typename T>
	size_t KDTreeBalancedPtInLeavesStackOpenCL<T>::getTreeDepth(size_t elCount) const
	{
		if (elCount <= 1)
			return 0;
		elCount --;
		size_t i = 31;
		for (; i >= 0; --i)
		{
			if (elCount & (1 << i))
				break;
		}
		return i+1;
	}

	template<typename T>
	void KDTreeBalancedPtInLeavesStackOpenCL<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos, const Vector minValues, const Vector maxValues)
	{
		const size_t count(last - first);
		//cerr << count << endl;
		if (count == 1)
		{
			const int dim = -2-(first->index);
			assert(pos < nodes.size());
			nodes[pos] = Node(dim);
			return;
		}
		
		// find the largest dimension of the box
		size_t cutDim = argMax<T>(maxValues - minValues);
		
		// compute number of elements
		const size_t rightCount(count/2);
		const size_t leftCount(count - rightCount);
		assert(last - rightCount == first + leftCount);
		
		// sort
		nth_element(first, first + leftCount, last, CompareDim(cutDim));
		
		// set node
		const T cutVal((first+leftCount)->pos.coeff(cutDim));
		nodes[pos] = Node(cutDim, cutVal);
		
		//cerr << pos << " cutting on " << cutDim << " at " << (first+leftCount)->pos[cutDim] << endl;
		
		// update bounds for left
		Vector leftMaxValues(maxValues);
		leftMaxValues[cutDim] = cutVal;
		// update bounds for right
		Vector rightMinValues(minValues);
		rightMinValues[cutDim] = cutVal;
		
		// recurse
		buildNodes(first, first + leftCount, childLeft(pos), minValues, leftMaxValues);
		buildNodes(first + leftCount, last, childRight(pos), rightMinValues, maxValues);
	}

	template<typename T>
	KDTreeBalancedPtInLeavesStackOpenCL<T>::KDTreeBalancedPtInLeavesStackOpenCL(const Matrix& cloud, const bool tryGPU):
	NearestNeighbourSearch<T>::NearestNeighbourSearch(cloud)
	{
		// build point vector and compute bounds
		BuildPoints buildPoints;
		buildPoints.reserve(cloud.cols());
		for (int i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			buildPoints.push_back(BuildPoint(v, i));
			const_cast<Vector&>(minBound) = minBound.cwise().min(v);
			const_cast<Vector&>(maxBound) = maxBound.cwise().max(v);
		}
		
		// create nodes
		nodes.resize(getTreeSize(cloud.cols()));
		buildNodes(buildPoints.begin(), buildPoints.end(), 0, minBound, maxBound);
		const unsigned maxStackDepth(getTreeDepth(nodes.size())*2 + 1);
		
		// looking for platforms, AMD drivers do not like the default for creating context
		vector<cl::Platform> platforms;
		cl::Platform::get(&platforms);
		for(vector<cl::Platform>::iterator i = platforms.begin(); i != platforms.end(); ++i)
		{
			cerr << "platform " << i - platforms.begin() << " is " << (*i).getInfo<CL_PLATFORM_VENDOR>() << endl;
		}
		cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0 };
		// TODO: add a way to specify the platform
		
		// create OpenCL contexts
		if (tryGPU)
		{
			try {
				context = cl::Context(CL_DEVICE_TYPE_GPU, cps);
			} catch (cl::Error e) {
				cerr << "Cannot find GPU for OpenCL" << endl;
				context = cl::Context(CL_DEVICE_TYPE_CPU, cps);
			}
		}
		else
			context = cl::Context(CL_DEVICE_TYPE_CPU, cps);
		
		// build and load source files
		cl::Program::Sources sources;
		// build defines
		ostringstream oss;
		oss << "typedef " << TypeName<T>::name << " T\n";
		oss << "#define DIM_COUNT " << cloud.rows() << "\n";
		oss << "#define POINT_STRIDE " << cloud.stride() << "\n";
		oss << "#define MAX_K " << MAX_K << "\n";
		oss << "#define MAX_STACK_DEPTH " << maxStackDepth << "\n";
		cerr << "params:\n" << oss.str() << endl;
		const size_t defLen(oss.str().length());
		char *defContent(new char[defLen+1]);
		strcpy(defContent, oss.str().c_str());
		sources.push_back(std::make_pair(defContent, defLen + 1));
		// load files
		const char* files[] = { "structure.cl", "knn.cl", NULL };
		for (const char** file = files; *file != NULL; ++file) {
			std::ifstream stream(*file);
			if (!stream.good())
				throw runtime_error((string("cannot open file: ") + *file));
			
			stream.seekg(0, std::ios_base::end);
			size_t size(stream.tellg());
			stream.seekg(0, std::ios_base::beg);
			
			char* content(new char[size + 1]);
			std::copy(std::istreambuf_iterator<char>(stream),
					  std::istreambuf_iterator<char>(), content);
			content[size] = '\0';
			
			sources.push_back(std::make_pair(content, size));
		}
		cl::Program program(context, sources);
		
		// build
		std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
		cl::Error error(CL_SUCCESS);
		try {
			program.build(devices);
		} catch (cl::Error e) {
			error = e;
		}
		
		// dump
		for (cl::Devices::const_iterator it = devices.begin(); it != devices.end(); ++it) {
			cerr << "device : " << it->getInfo<CL_DEVICE_NAME>() << "\n";
			cerr << "compilation log:\n" << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(*it) << endl;
		}
		// cleanup sources
		for (cl::Program::Sources::iterator it = sources.begin(); it != sources.end(); ++it) {
			delete[] it->first;
		}
		sources.clear();
		
		// make sure to stop if compilation failed
		if (error.err() != CL_SUCCESS)
			throw error;
		
		// build kernel and command queue
		knnKernel = cl::Kernel(program, "knnKDTree");
		queue = cl::CommandQueue(context, devices.back());
		
		// map nodes
		const size_t nodesCLSize(nodes.size() * sizeof(Node));
		nodesCL = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, nodesCLSize, &nodes);
		knnKernel.setArg(0, nodesCL);// TODO: manage alignment, see sect 6.1.5 
		// map cloud
		if (!(cloud.Flags & Eigen::DirectAccessBit) || (cloud.Flags & Eigen::RowMajorBit))
			throw runtime_error("wrong memory mapping in point cloud");
		const size_t cloudCLSize(cloud.cols() * cloud.stride());
		cloudCL = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, cloudCLSize, const_cast<T*>(&cloud.coeff(0,0)));
		knnKernel.setArg(1, cloudCL);
		
		// TODO: optimise by caching context and programs
	}
	
	
	template<typename T>
	typename KDTreeBalancedPtInLeavesStackOpenCL<T>::IndexMatrix KDTreeBalancedPtInLeavesStackOpenCL<T>::knnM(const Matrix& query, const Index k, const T epsilon, const unsigned optionFlags) 
	{
		// check K
		if (k > MAX_K)
			throw runtime_error("number of neighbors too large for OpenCL");
		
		// check consistency of query wrt cloud
		if (query.stride() != cloud.stride() ||
			query.rows() != cloud.rows())
			throw runtime_error("query is not of the same dimensionality as the point cloud");
		// map query
		if (!(query.Flags & Eigen::DirectAccessBit) || (query.Flags & Eigen::RowMajorBit))
			throw runtime_error("wrong memory mapping in query data");
		const size_t queryCLSize(query.cols() * query.stride());
		cl::Buffer queryCL(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, queryCLSize, const_cast<T*>(&query.coeff(0,0)));
		knnKernel.setArg(2, queryCL);
		// map result
		IndexMatrix result(k, query.cols());
		assert((query.Flags & Eigen::DirectAccessBit) && (!(query.Flags & Eigen::RowMajorBit)));
		const int indexStride(result.stride());
		const size_t resultCLSize(result.cols() * indexStride);
		cl::Buffer resultCL(context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, resultCLSize, &result.coeffRef(0,0));
		knnKernel.setArg(3, resultCL);
		
		// set resulting parameters
		knnKernel.setArg(4, k);
		knnKernel.setArg(5, 1 + epsilon);
		knnKernel.setArg(6, optionFlags);
		knnKernel.setArg(7, indexStride);
		
		// execute query
		queue.enqueueNDRangeKernel(knnKernel, cl::NullRange, cl::NDRange(nodes.size()), cl::NullRange);
		queue.enqueueMapBuffer(resultCL, true, CL_MAP_READ, 0, result.cols() * indexStride, 0, 0);
		queue.finish();
		
		// TODO: if requested, sort results
		return result;
	}

	template<typename T>
	typename KDTreeBalancedPtInLeavesStackOpenCL<T>::IndexVector KDTreeBalancedPtInLeavesStackOpenCL<T>::knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags)
	{
		Matrix m(query.size(), 1);
		m = query;
		return knnM(m, k, epsilon, optionFlags).col(0);
	}

	template struct KDTreeBalancedPtInLeavesStackOpenCL<float>;
	template struct KDTreeBalancedPtInLeavesStackOpenCL<double>;
	
	/*
	
	template<typename T>
	size_t KDTreeBalancedPtInLeavesStackOpenCL<T>::getTreeSize(size_t elCount) const
	{
		// FIXME: 64 bits safe stuff, only work for 2^32 elements right now
		assert(elCount > 0);
		elCount --;
		size_t count = 0;
		int i = 31;
		for (; i >= 0; --i)
		{
			if (elCount & (1 << i))
				break;
		}
		for (int j = 0; j <= i; ++j)
			count |= (1 << j);
		count <<= 1;
		count |= 1;
		return count;
	}

	template<typename T>
	void KDTreeBalancedPtInLeavesStackOpenCL<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos, const Vector minValues, const Vector maxValues)
	{
		const size_t count(last - first);
		//cerr << count << endl;
		if (count == 1)
		{
			const int dim = -2-(first->index);
			assert(pos < nodes.size());
			nodes[pos] = Node(dim);
			return;
		}
		
		// find the largest dimension of the box
		size_t cutDim = argMax<T>(maxValues - minValues);
		
		// compute number of elements
		const size_t rightCount(count/2);
		const size_t leftCount(count - rightCount);
		assert(last - rightCount == first + leftCount);
		
		// sort
		nth_element(first, first + leftCount, last, CompareDim(cutDim));
		
		// set node
		const T cutVal((first+leftCount)->pos.coeff(cutDim));
		nodes[pos] = Node(cutDim, cutVal);
		
		//cerr << pos << " cutting on " << cutDim << " at " << (first+leftCount)->pos[cutDim] << endl;
		
		// update bounds for left
		Vector leftMaxValues(maxValues);
		leftMaxValues[cutDim] = cutVal;
		// update bounds for right
		Vector rightMinValues(minValues);
		rightMinValues[cutDim] = cutVal;
		
		// recurse
		buildNodes(first, first + leftCount, childLeft(pos), minValues, leftMaxValues);
		buildNodes(first + leftCount, last, childRight(pos), rightMinValues, maxValues);
	}

	template<typename T>
	KDTreeBalancedPtInLeavesStackOpenCL<T>::KDTreeBalancedPtInLeavesStackOpenCL(const Matrix& cloud, const bool tryGPU):
	NearestNeighbourSearch<T>::NearestNeighbourSearch(cloud)
	{
		// build point vector and compute bounds
		BuildPoints buildPoints;
		buildPoints.reserve(cloud.cols());
		for (int i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			buildPoints.push_back(BuildPoint(v, i));
			const_cast<Vector&>(minBound) = minBound.cwise().min(v);
			const_cast<Vector&>(maxBound) = maxBound.cwise().max(v);
		}
		
		// create nodes
		nodes.resize(getTreeSize(cloud.cols()));
		buildNodes(buildPoints.begin(), buildPoints.end(), 0, minBound, maxBound);
		//for (size_t i = 0; i < nodes.size(); ++i)
		//	cout << i << ": " << nodes[i].dim << " " << nodes[i].cutVal << endl;
	}

	template<typename T>
	typename KDTreeBalancedPtInLeavesStackOpenCL<T>::IndexVector KDTreeBalancedPtInLeavesStackOpenCL<T>::knn(const Vector& query, const Index k, const T epsilon, const unsigned optionFlags)
	{
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T>::ALLOW_SELF_MATCH);
		
		assert(nodes.size() > 0);
		Heap heap(k);
		Vector off(Vector::Zero(query.size()));
		
		statistics.lastQueryVisitCount = 0;
		
		recurseKnn(query, 0, 0, heap, off, 1 + epsilon, allowSelfMatch);
		
		if (optionFlags & NearestNeighbourSearch<T>::SORT_RESULTS)
			heap.sort();
		
		statistics.totalVisitCount += statistics.lastQueryVisitCount;
		
		return heap.getIndexes();
	}

	template<typename T>
	void KDTreeBalancedPtInLeavesStackOpenCL<T>::recurseKnn(const Vector& query, const size_t n, T rd, Heap& heap, Vector& off, const T maxError, const bool allowSelfMatch)
	{
		const Node& node(nodes[n]);
		const int cd(node.dim);
		
		++statistics.lastQueryVisitCount;
		
		if (cd < 0)
		{
			if (cd == -1)
				return;
			const int index(-(cd + 2));
			const T dist(dist2<T>(query, cloud.col(index)));
			if ((dist < heap.headValue()) &&
				(allowSelfMatch || (dist > numeric_limits<T>::epsilon()))
			)
				heap.replaceHead(index, dist);
		}
		else
		{
			const T old_off(off.coeff(cd));
			const T new_off(query.coeff(cd) - node.cutVal);
			if (new_off > 0)
			{
				recurseKnn(query, childRight(n), rd, heap, off, maxError, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd * maxError < heap.headValue())
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childLeft(n), rd, heap, off, maxError, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
			else
			{
				recurseKnn(query, childLeft(n), rd, heap, off, maxError, allowSelfMatch);
				rd += - old_off*old_off + new_off*new_off;
				if (rd * maxError < heap.headValue())
				{
					off.coeffRef(cd) = new_off;
					recurseKnn(query, childRight(n), rd, heap, off, maxError, allowSelfMatch);
					off.coeffRef(cd) = old_off;
				}
			}
		}
	}

	template struct KDTreeBalancedPtInLeavesStackOpenCL<float>;
	template struct KDTreeBalancedPtInLeavesStackOpenCL<double>;
	*/
	
	//@}
}

#endif // HAVE_OPENCL
