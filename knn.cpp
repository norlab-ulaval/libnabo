#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

template<typename T>
T dist2(const typename Eigen::Matrix<T, Eigen::Dynamic, 1> v0, const typename Eigen::Matrix<T, Eigen::Dynamic, 1> v1)
{
	return (v0 - v1).squaredNorm();
}
	

template<typename T>
struct KDTree
{
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;

protected:
	struct BuildPoint
	{
		Vector pos;
		size_t index;
		BuildPoint(const Vector& pos =  Vector(), const size_t index = 0): pos(pos), index(index) {}
	};
	typedef vector<BuildPoint> BuildPoints;
	typedef typename BuildPoints::iterator BuildPointsIt;
	typedef typename BuildPoints::const_iterator BuildPointsCstIt;
	
	struct CompareDim
	{
		size_t dim;
		CompareDim(const size_t dim):dim(dim){}
		bool operator() (const BuildPoint& p0, const BuildPoint& p1) { return p0.pos(dim) < p1.pos(dim); }
	};
	
	struct SearchElement
	{
		size_t index;
		T minDist;
		
		SearchElement(const size_t index, const T minDist): index(index), minDist(minDist) {}
		friend bool operator<(const SearchElement& e0, const SearchElement& e1) { return e0.minDist < e1.minDist; }
	};
	
	// TODO: add index
	struct Node
	{
		Vector pos;
		int dim; // -1 == leaf, -2 == invalid
		Node(const Vector& pos = Vector(), const int dim = -2):pos(pos), dim(dim) {}
	};
	typedef vector<Node> Nodes;
	
	const size_t dim;
public:
	Vector minBound;
	Vector maxBound;
protected:
	Nodes nodes;
	
	inline size_t childLeft(size_t pos) const { return 2*pos + 1; }
	inline size_t childRight(size_t pos) const { return 2*pos + 2; }
	inline size_t parent(size_t pos) const { return (pos-1)/2; }
	size_t getTreeSize(size_t size) const;
	size_t argMax(const Vector& v) const;
	void buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos);
	void dump(const Vector minValues, const Vector maxValues, const size_t pos) const;
	
public:
	KDTree(const Matrix& cloud);
	Vector nn(const Vector& query);
};

template<typename T>
size_t KDTree<T>::getTreeSize(size_t elCount) const
{
	// FIXME: 64 bits safe stuff, only work for 2^32 elements right now
	size_t count = 0;
	int i = 31;
	for (; i >= 0; --i)
	{
		if (elCount & (1 << i))
			break;
	}
	for (int j = 0; j <= i; ++j)
		count |= (1 << j);
	//cerr << "tree size " << count << " (" << elCount << " elements)\n";
	return count;
}

template<typename T>
size_t KDTree<T>::argMax(const Vector& v) const
{
	T maxVal(0);
	size_t maxIdx;
	for (size_t i = 0; i < dim; ++i)
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
void KDTree<T>::buildNodes(const BuildPointsIt first, const BuildPointsIt last, const size_t pos)
{
	const size_t count(last - first);
	//cerr << count << endl;
	if (count == 1)
	{
		nodes[pos] = Node(first->pos, -1);
		return;
	}
	
	// estimate variance
	// get mean
	Vector mean(Vector::Zero(dim));
	for (BuildPointsCstIt it(first); it != last; ++it)
		mean += it->pos;
	mean /= last - first;
	// get sum of variance
	Vector var(Vector::Zero(dim));
	for (BuildPointsCstIt it(first); it != last; ++it)
		var += (it->pos - mean).cwise() * (it->pos - mean);
	// get dimension of maxmial variance
	const size_t cutDim = argMax(var);
	
	// sort
	sort(first, last, CompareDim(cutDim));
	
	// set node
	const size_t recurseCount(count-1);
	const size_t rightCount(recurseCount/2);
	const size_t leftCount(recurseCount-rightCount);
	assert(last - rightCount == first + leftCount + 1);
	
	nodes[pos] = Node((first+leftCount)->pos, cutDim);
	
	//cerr << pos << " cutting on " << cutDim << " at " << (first+leftCount)->pos[cutDim] << endl;
	
	// recurse
	if (count > 2)
	{
		buildNodes(first, first + leftCount, childLeft(pos));
		buildNodes(first + leftCount + 1, last, childRight(pos));
	}
	else
	{
		nodes[childLeft(pos)] = Node(first->pos, -1);
		nodes[childRight(pos)] = Node(Vector(), -2);
	}
}

template<typename T>
void KDTree<T>::dump(const Vector minValues, const Vector maxValues, const size_t pos) const
{
	const Node& node(nodes[pos]);
	
	if (node.dim >= -1)
	{
		if (dim == 2)
			cout << "<circle cx=\"" << 100*node.pos(0) << "\" cy=\"" << 100*node.pos(1) << "\" r=\"1\" stroke=\"black\" stroke-width=\"0.2\" fill=\"red\"/>" << endl;
		else
			cout << "pt at\n" << node.pos << endl;
	}
	if (node.dim >= 0)
	{
		//cerr << "in bounds:\n" << minValues << "\nto\n" << maxValues << endl;
		
		// update bounds for left
		Vector leftMaxValues(maxValues);
		leftMaxValues[node.dim] = node.pos[node.dim];
		// update bounds for right
		Vector rightMinValues(minValues);
		rightMinValues[node.dim] = node.pos[node.dim];
		
		// print line
		if (dim == 2)
			cout << "<line x1=\"" << 100*rightMinValues(0) << "\" y1=\"" << 100*rightMinValues(1) << "\" x2=\"" << 100*leftMaxValues(0) << "\" y2=\"" << 100*leftMaxValues(1) << "\" style=\"stroke:rgb(0,0,0);stroke-width:0.2\"/>" << endl;
		else
			cout << "cut from\n" << rightMinValues << "\nto\n" << leftMaxValues << endl;
		// recurs
		dump(minValues, leftMaxValues, childLeft(pos));
		dump(rightMinValues, maxValues, childRight(pos));
	}
}

template<typename T>
KDTree<T>::KDTree(const typename KDTree<T>::Matrix& cloud):
	dim(cloud.rows()),
	minBound(Vector::Constant(dim, numeric_limits<T>::max())),
	maxBound(Vector::Constant(dim, numeric_limits<T>::min()))
{
	// build point vector and compute bounds
	vector<BuildPoint> buildPoints;
	for (size_t i = 0; i < cloud.cols(); ++i)
	{
		const Vector& v(cloud.col(i));
		buildPoints.push_back(BuildPoint(v, i));
		minBound = minBound.cwise().min(v);
		maxBound = maxBound.cwise().max(v);
	}
	
	// create nodes
	nodes.resize(getTreeSize(cloud.cols()));
	buildNodes(buildPoints.begin(), buildPoints.end(), 0);
	
	// dump nodes
	//dump(minBound, maxBound, 0);
}

template<typename T>
typename KDTree<T>::Vector KDTree<T>::nn(const Vector& query)
{
	typedef priority_queue<SearchElement> Queue;
	Queue queue;
	
	queue.push(SearchElement(0, 0));
	
	T bestDist(numeric_limits<T>::max()); 
	size_t bestIndex(0);
	size_t visitCount(0);
	while (!queue.empty())
	{
		SearchElement el(queue.top());
		queue.pop();
		
		// nothing is closer, we found best
		if (el.minDist > bestDist)
			break;
		
		const Node& node(nodes[el.index]);
		// TODO: optimise, do not push these in first place
		if (node.dim == -2)
			continue;
		
		const T dist(dist2(node.pos, query));
		if (dist < bestDist)
		{
			bestIndex = el.index;
			bestDist = dist;
		}
		
		if (node.dim < 0)
			continue;
		
		const T offset(query[node.dim] - node.pos[node.dim]);
		const T offset2(offset * offset);
		if (offset > 0)
		{
			if (offset2 < bestDist)
				queue.push(SearchElement(childLeft(el.index), offset2));
			queue.push(SearchElement(childRight(el.index), 0));
		}
		else
		{
			queue.push(SearchElement(childLeft(el.index), 0));
			if (offset2 < bestDist)
				queue.push(SearchElement(childRight(el.index), offset2));
		}
		// TODO: optimise using loop there to relieve queue use
		
		++visitCount;
	}
	cerr << "visit count kdtree: " << visitCount << endl;
	return nodes[bestIndex].pos;
}



// brute force nearest neighbor for comparison
template<typename T>
typename KDTree<T>::Vector bfnn(const typename KDTree<T>::Vector& query, const typename KDTree<T>::Matrix data)
{
	T bestDist(numeric_limits<T>::max()); 
	size_t bestIndex(0);
	
	for (size_t i = 0; i < data.cols(); ++i)
	{
		const T dist(dist2<float>(data.col(i), query));
		if (dist < bestDist)
		{
			bestIndex = i;
			bestDist = dist;
		}
	}
	
	cerr << "visit count bf: " << data.cols() << endl;
	
	return data.col(bestIndex);
}

template<typename T>
typename KDTree<T>::Matrix load(const char *fileName)
{
	typedef typename KDTree<T>::Matrix Matrix;
	
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



int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cerr << "Usage " << argv[0] << " DATA" << endl;
		return 1;
	}
	
	/*cout << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
	cout << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;

	cout << "<svg width=\"100%\" height=\"100%\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl;
	*/
	
	KDTree<float>::Matrix d(load<float>(argv[1]));
	//cerr << d.rows() << " " << d.cols() << endl;
	KDTree<float> kdtree(d);
	typedef KDTree<float>::Vector Vector;
	
	for (int i = 0; i < 10; ++i)
	{
		Vector q(kdtree.minBound + Vector::Random(kdtree.minBound.size()).cwise() * (kdtree.maxBound - kdtree.minBound));
		Vector v_bf(bfnn<float>(q, d));
		Vector v_kdtree(kdtree.nn(q));
		cout << "query:\n" << q << "\nbf:\n" << v_bf << "\nkdtree:\n" << v_kdtree << "\n\n";
		assert(dist2(v_bf, v_kdtree) < numeric_limits<float>::epsilon());
	}
		
	//cout << "</svg>" << endl;
	
	return 0;
}