/**
 * Defines the basic methods for a shift invariant space
 */
#include <sisl/primitives.hpp>
#include <functional>
#include <vector>

#ifndef __LATTICE__
#define __LATTICE__

namespace sisl {
template<class L, class GF, class O = float, class I = float> class filter;

template <class GF, class O = float, class I = float>
class shift_invariant_space3 {
public:
	shift_invariant_space3(){} // should initialize with h = 1
	shift_invariant_space3(const I &h){}
	~shift_invariant_space3(){}

	virtual void SV(const int &x, const int &y, const int &z, const O &value) = 0;
	virtual O GV(const int &x, const int &y, const int &z) const = 0;

	virtual void SN(const int &x, const int &y, const int &z, const vector3<O> &value) = 0;
	virtual vector3<O> GN(const int &x, const int &y, const int &z) = 0;

	virtual const int lIndex(const int &x, const int &y, const int &z) const = 0;
	virtual const int numberOfLatticeSites() const = 0;

	virtual I getScale() const = 0;
	virtual const int getResolution() const = 0;
	virtual vector3<I> getSitePosition(const int &x, const int &y, const int &z) const = 0;
	virtual vector3<int> getNearestIndex(const vector3<I> &) const = 0;

	// Evaluates the function  
	virtual const O f(const I &x, const I &y, const I &z) = 0;
	virtual vector3<O> grad_f(const I &x, const I &y, const I &z) = 0; 
	virtual const O f(const vector3<I>&) = 0;
	virtual vector3<O> grad_f(const vector3<I>&) = 0; 

	//
	virtual void forEachLatticeSite(std::function<O(const int &, const int &, const int &)>) = 0; //
	virtual void spatialFilter(std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int &)>)>) = 0; //
	virtual void frequencyFilter(std::function<O(const I &u, const I &v, const I &w, const I &h)> lambda) = 0; //

	//
	virtual std::string getLatticeName() const = 0;
	virtual std::string getBasisName() const = 0;

	// Convolution of two functions
	//inline shift_invariant_space3<GF, O, I> operator*(const shift_invariant_space3<GF, O, I> &v) {throw;};

	// Copy a function
	//inline shift_invariant_space3<GF, O, I> operator=(const shift_invariant_space3<GF, O, I> &v) {throw;};
};

template <class O = float, class I = float>
static std::function<O(const I &u, const I &v, const I &w, const I &h)> 
combineFilters(std::function<O(const I &u, const I &v, const I &w, const I &h)> f1, std::function<O(const I &u, const I &v, const I &w, const I &h)> f2) {
	return [&](const I &u, const I &v, const I &w, const I &h) ->  O {
		return f1(u,v,w,h)*f2(u,v,w,h);
	};
}

template <class GF, class O = float, class I = float>
class shift_invariant_space2 {
public:
	shift_invariant_space2(){} // should initialize with h = 1
	shift_invariant_space2(const I &h){}
	~shift_invariant_space2(){}

	virtual void SV(const int &x, const int &y, const int &z, const O &value) = 0;
	virtual O GV(const int &x, const int &y, const int &z) const = 0;
	virtual const int lIndex(const int &x, const int &y, const int &z) const = 0;

	virtual I getScale() const = 0;
	virtual vector3<I> getSitePosition(const int &x, const int &y, const int &z) const = 0;

	// Evaluates the function  
	virtual const O f(const I &x, const I &y, const I &z) = 0;
	virtual vector3<I> grad_f(const I &x, const I &y, const I &z) = 0; 
	virtual const O f(const vector3<I>&) = 0;
	virtual vector3<I> grad_f(const vector3<I>&) = 0; 

	//
	virtual void spatialFilter(O (*lambda)(const int &, const int &, O (*a)(const int &, const int &))) = 0; //
	virtual void frequencyFilter(O (*lamda)(const I &, const I &)) = 0; //

	// Convolution of two functions
	//inline shift_invariant_space2<GF, O, I> operator*(const shift_invariant_space2<GF, O, I> &v) {throw;};

	// Copy a function
	//inline shift_invariant_space2<GF, O, I> operator=(const shift_invariant_space2<GF, O, I> &v) {throw;};
};
};

#define MODP(x, d) (((x)%(d) + (d)) % (d))

#endif // __LATTICE__
