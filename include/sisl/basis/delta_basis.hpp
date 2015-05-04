#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _DELTA_KRONEKER_H_
#define _DELTA_KRONEKER_H_

namespace sisl{
template <class O, class I>
class delta_basis: public basis_function <O,I> {
private:
	static inline O bspline3(const I &t){
		if(t <= -2  || t > 2)
			return 0;
		if(t <= -1 && t > -2)
			return (4./3.) + t*(2. + t*(1. + (1./6.)*t));
		if(t <= 0 && t > -1)
			return (2./3.) + t*t*(-1. - 0.5*t);
		if( t <= 1 && t > 0)
			return (2./3.) + t*t*(-1. + 0.5*t);
		return (4./3.) + t*(-2. + t*(1. - (1./6.)*t));
	}

public:
	static std::string getBasisName(){
		return std::string("delta_basis");
	}
	static const O M(const I &x, const I &y, const I &z) { 
		if(x == y && y == z && z== x && z == 0.) 
			return 1.; 
		return 0; 
	}

	static const O M(const vector3<I> &p) { return (O)bspline3(p.i)*bspline3(p.j)*bspline3(p.k); 
		if(p.i == p.j && p.j == p.k && p.k == p.i && p.i == 0.) 
			return 1.; 
		return 0; 
	}

	// This function should return the intersection of the closure of the support of
	// the generator and the lattice.
	static std::vector<std::tuple<int,int,int>> getSupport() {
		using namespace std;
		std::vector<std::tuple<int,int,int>> support;
		support.push_back(make_tuple(0, 0, 0));
		return support;
	};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int>> getEffectiveSupport(const vector3<I> &p) {
		using namespace std;
		std::vector<std::tuple<int,int,int>> support;
		support.push_back(make_tuple(0, 0, 0));
		return support;
	}

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi2Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> autoCorrelation(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};

	inline static const O convolutionSum(const vector3<I> &p, const shift_invariant_space3< delta_basis, O, I> *lattice) {
				O sum = 0;

		vector3<I> vox;
		
		I dh = lattice->getScale();

		vox.i = p.i/dh;
		vox.j = p.j/dh;
		vox.k = p.k/dh;

		int 
			vx = (int)floor(vox.i), 
			vy = (int)floor(vox.j), 
			vz = (int)floor(vox.k);

		return lattice->GV(vx, vy, vz);
	}

	inline static const vector3<O> convolutionSumNormal(const vector3<I> &p, shift_invariant_space3< delta_basis, O, I> *lattice) {

		vector3<O> sum = vector3<O>(0,0,0);
		vector3<I> vox;
		
		I dh = lattice->getScale();

		vox.i = p.i/dh;
		vox.j = p.j/dh;
		vox.k = p.k/dh;

		int 
			vx = (int)floor(vox.i), 
			vy = (int)floor(vox.j), 
			vz = (int)floor(vox.k);

		return lattice->GN(vx, vy, vz);
	}
};
};
#endif //_TP_BSPLINE3_H_