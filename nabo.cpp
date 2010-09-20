#include "nabo.h"
#include <limits>
#include <algorithm>

namespace Nabo
{
	using namespace std;
	
	template<typename T>
	NearestNeighborSearch<T>::NearestNeighborSearch(const Matrix& cloud):
		cloud(cloud),
		dim(cloud.rows()),
		minBound(Vector::Constant(dim, numeric_limits<T>::max())),
		maxBound(Vector::Constant(dim, numeric_limits<T>::min()))
	{
		/*// build point vector and compute bounds
		vector<BuildPoint> buildPoints;
		for (size_t i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.col(i));
			buildPoints.push_back(BuildPoint(v, i));
			minBound = minBound.cwise().min(v);
			maxBound = maxBound.cwise().max(v);
		}*/
	}
	
	template struct NearestNeighborSearch<float>;
	template struct NearestNeighborSearch<double>;
}
