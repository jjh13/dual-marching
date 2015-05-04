/*
	1/2^n hyperspace.
*/

#ifndef __NVAR_HALFSPACE__
#define __NVAR_HALFSPACE__

#include <vector>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>

namespace sisl {
	template<int N, class I=double>
	class halfspace{
	public:
		halfspace(const Eigen::Matrix<I, Dynamic, Dynamic> &M) {
			this->T = M;
			this->Tinv = M.inverse();
		}

		bool in_halfspace(std::vector<I> x) {
			Eigen::Matrix<I, Dynamic, 1> v(x.size());
			for(int i=0; i < x.size(); i++) v(i,0) = x[i];

			v = Tinv * v; 
			for(int i=0; i < x.size(); i++) 
				if(v(i) < (I)0)
					return false;
			return true;

		}
	private:
		Eigen::Matrix<I, Dynamic, Dynamic> T, Tinv;
	};
}

#endif //__NVAR_HALFSPACE__