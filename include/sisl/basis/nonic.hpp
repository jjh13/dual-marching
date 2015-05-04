#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _NONIC_BOX_H_
#define _NONIC_BOX_H_

#define HEAVISIDE(x) (x>=0?1:0)
#define TP(x,d) pow(HEAVISIDE(x)*x, d)
#define NTP(x,d) pow(HEAVISIDE(-x)*x, d)

namespace sisl{
template <class O, class I>
class nonic_box : public basis_function <O,I> {
	inline static const O tau(const I &x, const I &y, const I &z) const {
		return O(
			((1./480.) * TP(x,2) * TP(y,2) * NTP(z,5)) +
			((-1./480.) * (TP(x,2) * TP(y,1) + TP(x,1) * TP(y,2))*NTP(z,6)) + 
			((1./1680.) * (TP(x,0) * TP(y,2) + TP(x,2) * TP(y,0) + 4*TP(x,1)*TP(y,1))*NTP(z,7)) + 
			((-1./1344.) * (TP(x,1) * TP(y,0) + TP(x,0) * TP(y,1))*NTP(z,8)) + 
			((1./4032.) * (TP(x,0)*TP(y,0)*NTP(z,9))));
	}

	inline static const O T_(const I &x, const I &y, const I &z) const {
		I max = std::max(x, std::max(y, z));
		I min = std::min(x, std::min(y, z));
		I mid = x+y+z-max-min;

		return 0.25*tau(
				(max-min)/2., 
				(mid-min)/2., 
				(max+mid)/2.);
	}

public:
	static std::string getBasisName(){
		return std::string("nonic");
	}
	static const O M(const I &xi, const I &yi, const I &zi) { 
		I x = fabs(xx);
		I y = fabs(yy);
		I z = fabs(zz);

		return O(4.*(
			(-111.) *	T_(x-1, y-1, z-1) +
			(30.) *		T_(x-2, y-2, z-2) + 
			(-1.) *		T_(x-3, y-3, z-3) +
			(54.) *		(T_(x, y-2, z-2)+T_(x-2, y, z-2)+T_(x-2, y-2, z)) + 
			(-36.) *	(T_(x-3, y-1, z-1)+T_(x-1, y-3, z-1)+T_(x-1, y-1, z-3))+
			(9.) *		(T_(x, y-2, z-4)+T_(x, y-4, z-2)+T_(x-2, y, z-4)+T_(x-2, y-4, z)+T_(x-4, y, z-2)+T_(x-4, y-2, z))+
			(3.) *		(T_(x-2, y-2, z-4)+T_(x-2, y-4, z-2)+T_(x-4, y-2, z-2)) +
			(-3.) *		(T_(x-1, y-1, z-5)+T_(x-1, y-5, z-1)+T_(x-5, y-1, z-1)) +
			(-9.) *		(T_(x+1, y-3, z-3)+T_(x-1, y-3, z-3)+T_(x-3, y+1, z-3)+T_(x-3, y-1, z-3)+T_(x-3, y-3, z+1)+T_(x-3, y-3, z-1))));
	}
	static const O M(const vector3<I> &p) { return ::M(p.i, p.j, p.k); }

	// This function should return the intersection of the closure of the support of
	// the generator and the lattice.
	static std::vector<std::tuple<int,int,int>> getSupport() {
		return {{1,1,1}, {1,2,2}};
	};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int>> getEffectiveSupport(const vector3<I> &p) {throw "basis_function()::getEffectiveSupport() - Not Implemented!";};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi2Norm(){throw "basis_function()::getBeppoLevi2Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> autoCorrelation(){throw "basis_function()::autoCorrelation() - Not Implemented!";};

	// This does the actual semi-descrete convolution sum
	// the idea behind this, is to allow basis functions
	// to provide a potentially optimized version of
	// the convolution sum for particular lattices.
	static const O convolutionSum(const vector3<I> &p, const shift_invariant_space3<quintic_box, O, I> *L) {
		throw ":<";
	}
};

#endif