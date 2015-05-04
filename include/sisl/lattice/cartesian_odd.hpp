/**
 * General purpose class for cartesian lattices with resolution 
 * i+2*j+2*k+2, and zero boundary conditions, on the unit cube.
 *
 * @author Joshua Horacsek
 *
 **/

#ifndef _CARTESIAN_H
#define _CARTESIAN_H

#include <sisl/array.hpp>
#include <sisl/sparse_array.hpp>
#include <sisl/lattice.hpp>
#include <sisl/basis_function.hpp>
#include <sisl/fftwalloc.hpp>
#include <fftw3.h>

#include <iostream>

namespace sisl{
template <class GF, class O = float, class I = float>
class cartesian_odd : public shift_invariant_space3<GF,O,I> {
public:
	cartesian_odd() : cartesian_odd(1) {} // should initialize with h = 1
	cartesian_odd(const I &h) : 
		dh(h), 
		res((1./h)+1),
		NA(res-2,res-2,res-2,vector3<O>(0,0,0)),
		L(res-2,res-2,res-2) {};
	~cartesian_odd(){}

	virtual O GV(const int &x, const int &y, const int &z) const {
		int xr = MODP(x, res-1);
		int yr = MODP(y, res-1);
		int zr = MODP(z, res-1);

		if(xr == 0 || yr == 0 || zr == 0) 
			return 0;

		return L(xr-1, yr-1, zr-1) *
					(x < 0 || x > (res - 1) ? -1 : 1) * 
					(y < 0 || y > (res - 1) ? -1 : 1) * 
					(z < 0 || z > (res - 1) ? -1 : 1);
	}

	virtual void SV(const int &x, const int &y, const int &z, const O &value) {
		int xr = MODP(x, (res - 1));
		int yr = MODP(y, (res - 1));
		int zr = MODP(z, (res - 1));

		if(xr == 0 || yr == 0 || zr == 0) 
			return;

		L(xr-1, yr-1, zr-1) = value;
	}
	virtual vector3<O> GN(const int &x, const int &y, const int &z) {
		int xr = MODP(x, res-1);
		int yr = MODP(y, res-1);
		int zr = MODP(z, res-1);

		if(xr == 0 || yr == 0 || zr == 0) 
			return vector3<O>(0,0,0);

		if(!NA.hasValue(xr-1,yr-1,zr-1))
			NA(xr-1,yr-1,zr-1) = GF::approximateGradient(x,y,z, this);

		return NA(xr-1, yr-1, zr-1) *
					(x < 0 || x > (res - 1) ? -1 : 1) * 
					(y < 0 || y > (res - 1) ? -1 : 1) * 
					(z < 0 || z > (res - 1) ? -1 : 1);
	}

	virtual void SN(const int &x, const int &y, const int &z, const vector3<O> &value) {
		int xr = MODP(x, (res - 1));
		int yr = MODP(y, (res - 1));
		int zr = MODP(z, (res - 1));

		if(xr == 0 || yr == 0 || zr == 0) 
			return;

		NA(xr-1, yr-1, zr-1) = value;
	}
	virtual const int lIndex(const int &x, const int &y, const int &z) const { 
		int xr = MODP(x, (res - 1));
		int yr = MODP(y, (res - 1));
		int zr = MODP(z, (res - 1));

		if(xr == 0 || yr == 0 || zr == 0 ||x >=  (res - 1) || y >=  (res - 1) || z >= (res - 1) || x < 0 || y < 0 || z < 0 ) 
			return -1;

		return L.lIndex(xr-1, yr-1, zr-1);
	}
	virtual const int numberOfLatticeSites() const {
		return L.size();
	}

	virtual I getScale() const {return dh;}
	virtual const int getResolution() const {
		return res-2;
	}

	virtual vector3<I> getSitePosition(const int &x, const int &y, const int &z) const { 
		return vector3<I>(I(x)*dh,I(y)*dh,I(z)*dh); 
	}
	virtual vector3<int> getNearestIndex(const vector3<I> &pos) const {
		return vector3<int>((int)floor(pos.i/dh), (int)floor(pos.j/dh), (int)floor(pos.k/dh));
	}

	// Evaluates the function  
	virtual const O f(const I &x, const I &y, const I &z) { return GF::convolutionSum(vector3<I>(x,y,z), this); }
	virtual vector3<O> grad_f(const I &x, const I &y, const I &z) { return GF::convolutionSumNormal(vector3<I>(x,y,z), this); }

	// Evaluates the function  
	virtual const O f(const vector3<I>& p ) {return GF::convolutionSum(p, this); }
	virtual vector3<O> grad_f(const vector3<I>& p) { return GF::convolutionSumNormal(p, this); }

	//
	virtual void forEachLatticeSite(std::function<O(const int &, const int &, const int &)> lambda) {
		#pragma omp parallel for
		for(int i = 1; i < this->res - 1; i++)
			for(int j = 1; j < this->res - 1; j++)
				for(int k = 1; k < this->res - 1; k++) {
					SV(i,j,k, lambda(i,j,k));
				}
	}
	
	//
	virtual void spatialFilter(std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> lambda) {
		cartesian_odd <GF, O, I> src = *this; // Copy the lattice

		std::function<O(const int &, const int &, const int &)> accessFunc = 
			[&](const int &x, const int &y, const int &z) ->  O {
				return (O) src.GV(x,y,z); 
			};
		
		#pragma omp parallel for
		for(int i = 1; i < this->res - 1; i++)
			for(int j = 1; i < this->res - 1; j++)
				for(int k = 1; i < this->res - 1; k++) {
					SV(i,j,k, lambda(i,j,k, accessFunc));
				}
	} 
	//
	virtual void frequencyFilter(std::function<O(const I &u, const I &v, const I &w, const I &h)> lambda) {
		I scale =  (0.125*dh*dh*dh);
		I invdh = 1./(dh*dh);
		O *lattice = L.getArray();

		// Do the forward transform
		fftw_plan p = fftw_plan_r2r_3d(res-2, res-2, res-2,
					lattice, lattice, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);

		// Apply the filter(s)
		#pragma omp parallel for
		for(int i = 0; i<res-2; i++) 
			for(int j = 0; j<res-2; j++)
				for(int k = 0; k<res-2; k++){
					I u = M_PI*I(i+1)*dh;
					I v = M_PI*I(j+1)*dh;
					I w = M_PI*I(k+1)*dh;
					L(i,j,k) *= lambda(u,v,w,dh) * scale;
				}
		
		// Do the inverse transform
		p = fftw_plan_r2r_3d(res-2, res-2, res-2,
					lattice, lattice, FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);

		fftw_execute(p);
		fftw_destroy_plan(p);
	} 

	// Convolution of two functions
	//inline cartesian_odd<GF, O, I> operator*(const cartesian_odd<GF, O, I> &v) {throw;};

	// Copy a function
	inline cartesian_odd<GF, O, I> operator=(const cartesian_odd<GF, O, I> &v) {throw;};

	virtual std::string getLatticeName() const {
		return std::string("cartesian");
	}
	virtual std::string getBasisName() const {
		return GF::getBasisName();
	}
private:
	I dh;
	int res;
	array3<O, _fftwalloc<O> > L;
	sparse_array3<vector3<O> > NA;
};
};

#endif // _CARTESIAN_H
