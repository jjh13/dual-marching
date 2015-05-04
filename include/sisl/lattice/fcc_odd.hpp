/**
 * General purpose class for cartesian lattices with resolution 
 * i+2*j+2*k+2, and zero boundary conditions, on the unit cube.
 *
 * @author Joshua Horacsek
 *
 **/

#ifndef _FCC_H
#define _FCC_H

#include <sisl/array.hpp>
#include <sisl/sparse_array.hpp>
#include <sisl/lattice.hpp>
#include <sisl/basis_function.hpp>
#include <sisl/fftwalloc.hpp>
#include <fftw3.h>

#include <iostream>

namespace sisl{
template <class GF, class O = float, class I = float>
class fcc_odd : public shift_invariant_space3<GF,O,I> {
public:
	fcc_odd() : fcc_odd(1) {} // should initialize with h = 1
	fcc_odd(const I &h) : 
		dh(h), 
		res((1./h)+1),
		NA(res,res,res/2, 4,vector3<O>(0,0,0)),
		L(res,res,res/2, 4) {};
	~fcc_odd(){}

	virtual O GV(const int &x, const int &y, const int &z) const  {
		int index =  ((x - 1) % 2) | (((y - 1) % 2) << 1) | (((z - 1) % 2) << 2);
		const char lookup[10] = { 0, -1, -1, -1, -1, 1, 2, -1, -1, 3, };
		if(lookup[index] > 3 || lookup[index] < 0) throw;
		return L((x-1) >> 1, (y-1) >> 1, (z-1) >> 1, lookup[index]);
	};
	virtual void SV(const int &x, const int &y, const int &z, const O &value)  {
		int index =  ((x - 1) % 2) | (((y - 1) % 2) << 1) | (((z - 1) % 2) << 2);
		const char lookup[10] = { 0, -1, -1, -1, -1, 1, 2, -1, -1, 3, };
		if(lookup[index] > 3 || lookup[index] < 0) throw;
		L((x-1) >> 1, (y-1) >> 1, (z-1) >> 1, lookup[index]) = value;
	};
	virtual vector3<O> GN(const int &x, const int &y, const int &z)  {
		int index =  ((x - 1) % 2) | (((y - 1) % 2) << 1) | (((z - 1) % 2) << 2);
		const char lookup[10] = { 0, -1, -1, -1, -1, 1, 2, -1, -1, 3, };
		if(lookup[index] > 3 || lookup[index] < 0) throw;
		return NA((x-1) >> 1, (y-1) >> 1, (z-1) >> 1, lookup[index]);
	};
	virtual void SN(const int &x, const int &y, const int &z, const vector3<O> &value)  {
		int index =  ((x - 1) % 2) | (((y - 1) % 2) << 1) | (((z - 1) % 2) << 2);
		const char lookup[10] = { 0, -1, -1, -1, -1, 1, 2, -1, -1, 3, };
		if(lookup[index] > 3 || lookup[index] < 0) throw;
		NA((x-1) >> 1, (y-1) >> 1, (z-1) >> 1, lookup[index]) = value;
	};
	virtual const int lIndex(const int &x, const int &y, const int &z) const  {
		int index =  ((x - 1) % 2) | (((y - 1) % 2) << 1) | (((z - 1) % 2) << 2);
		const char lookup[10] = { 0, -1, -1, -1, -1, 1, 2, -1, -1, 3, };
		if(lookup[index] > 3 || lookup[index] < 0) throw;
		return L.lIndex((x-1) >> 1, (y-1) >> 1, (z-1) >> 1, lookup[index]);
	};

	virtual const int numberOfLatticeSites() const {
		return L.size();
	}

	virtual I getScale() const  { return dh; };
	virtual const int getResolution() const  {return res;};

	virtual vector3<I> getSitePosition(const int &x, const int &y, const int &z) const  {
			// int xx = x*2 + ((y&1) ^ (z&1));
			return vector3<I>(I(x)*dh,I(y)*dh,I(z)*dh);

	}
	virtual vector3<int> getNearestIndex(const vector3<I> &pos) const  {throw;};

	// Evaluates the function  
	virtual const O f(const I &x, const I &y, const I &z) { return GF::convolutionSum(vector3<I>(x,y,z), this); }
	virtual vector3<O> grad_f(const I &x, const I &y, const I &z) { return GF::convolutionSumNormal(vector3<I>(x,y,z), this); }

	// Evaluates the function  
	virtual const O f(const vector3<I>& p ) {return GF::convolutionSum(p, this); }
	virtual vector3<O> grad_f(const vector3<I>& p) { return GF::convolutionSumNormal(p, this); }

	//
	virtual void forEachLatticeSite(std::function<O(const int &, const int &, const int &)> lambda) {
		#pragma omp parallel for
		for(unsigned int i = 2; i < res-2; i += 2)
			for(unsigned int j = 1; j < res - 1; j++)
				for(unsigned int k = 1; k < res - 1; k++) {

					int ii = i*2 + ((j&1) ^ (k&1));

					/* The little bit-hack operations encode the following
					if(j % 2) {
						if(k % 2){
							ii = i*2;
						}else{
							ii = i*2 + 1;
						}
					} else {
						if(k % 2){
							ii = i*2 + 1;
						}else{
							ii = i*2;
						}
					}*/
					SV(ii, j, k, lambda(ii, j, k));
				}

	}
	
	//
	virtual void spatialFilter(std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> lambda) {
		throw;
	} 
	//
	virtual void frequencyFilter(std::function<O(const I &u, const I &v, const I &w, const I &h)> lambda) {
		throw;
	} 

	// Convolution of two functions
	//inline fcc_odd<GF, O, I> operator*(const cartesian_odd<GF, O, I> &v) {throw;};

	// Copy a function
	inline fcc_odd<GF, O, I> operator=(const fcc_odd<GF, O, I> &v) {throw;};

	virtual std::string getLatticeName() const {
		return std::string("fcc");
	}
	virtual std::string getBasisName() const {
		return GF::getBasisName();
	}
private:
	I dh;
	int res;
	array4<O, _fftwalloc<O> > L;
	sparse_array4<vector3<O> > NA;
};
};

#endif // _CARTESIAN_H