
#include <sisl/array.hpp>
#include <sisl/sparse_array.hpp>
#include <sisl/lattice.hpp>
#include <sisl/basis_function.hpp>
#include <sisl/fftwalloc.hpp>
#include <fftw3.h>

#include <iostream>

#ifndef _BCC_ODD_H
#define _BCC_ODD_H

namespace sisl{
template <class GF, class O = float, class I = float>
class bcc_odd : public shift_invariant_space3<GF,O,I> {
public:
	bcc_odd() : bcc_odd(.5) {} // should initialize with h = 1
	bcc_odd(const I &h) : 
		dh(h), 
		res(1./(2.*dh)), 
		NA(res,res,res*2,vector3<O>(0,0,0)), 
		L(res,res,res*2){};
	~bcc_odd(){}

	virtual O GV(const int &x, const int &y, const int &z) const {
		int xr = (MODP(x, res * 2) - abs(z)%2)/2;
		int yr = (MODP(y, res * 2) - abs(z)%2)/2;
		int zr = MODP(z, res * 2);

		if( x % (res * 2) == 0 || y % (res * 2) == 0 || z == 0)
			return 0;

		return L(xr, yr, zr) *
			(x < 0 || x > (2*res) ? -1 : 1) *
			(y < 0 || y > (2*res) ? -1 : 1) *
			(z < 0 || z > (2*res) ? -1 : 1);
	}

	virtual void SV(const int &x, const int &y, const int &z, const O &value) {
		int xr = (MODP(x, res * 2) - abs(z)%2)/2;
		int yr = (MODP(y, res * 2) - abs(z)%2)/2;
		int zr = MODP(z, res * 2);

		L(xr, yr, zr) = value;
	}

	virtual vector3<O> GN(const int &x, const int &y, const int &z) {
		int xr = (MODP(x, res * 2) - abs(z)%2)/2;
		int yr = (MODP(y, res * 2) - abs(z)%2)/2;
		int zr = MODP(z, res * 2);

		if( x % (res * 2) == 0 || y % (res * 2) == 0 || z == 0)
			return vector3<O>(0,0,0);

		if(!NA.hasValue(xr, yr, zr))
			NA(xr, yr, zr) = GF::approximateGradient(x,y,z, this);
		
		return NA(xr, yr, zr) *
			(x < 0 || x > (2*res) ? -1 : 1) *
			(y < 0 || y > (2*res) ? -1 : 1) *
			(z < 0 || z > (2*res) ? -1 : 1);
	}

	virtual void SN(const int &x, const int &y, const int &z, const vector3<O> &value) {
		int xr = (MODP(x, res * 2) - abs(z)%2)/2;
		int yr = (MODP(y, res * 2) - abs(z)%2)/2;
		int zr = MODP(z, res * 2);

		NA(xr, yr, zr) = value;
	}
	virtual const int lIndex(const int &x, const int &y, const int &z) const { 
		int xr = (MODP(x, res * 2) - abs(z)%2)/2;
		int yr = (MODP(y, res * 2) - abs(z)%2)/2;
		int zr = MODP(z, res * 2);

		return L.lIndex(xr, yr, zr);
	}

	virtual const int numberOfLatticeSites() const {return L.size(); }
	virtual I getScale() const {return dh;}
	virtual const int getResolution() const {return res;}
	virtual vector3<I> getSitePosition(const int &x, const int &y, const int &z) const { return vector3<I>(I(x)*dh,I(y)*dh,I(z)*dh); }

	virtual vector3<int> getNearestIndex(const vector3<I> &pt) const {
			int P1[3],P2[3],P3[3],P4[3];
			int I1[3],I2[3];

			I x = pt.i/dh;
			I y = pt.j/dh;
			I z = pt.k/dh;

			sisl::vector3<I> BCCvox(
				(x + y) / 2,
				(x + z) / 2,
				(y + z) / 2);

			int vx = (int)floor(BCCvox.i),
				vy = (int)floor(BCCvox.j),
				vz = (int)floor(BCCvox.k);

			sisl::vector3<I> ga(
				BCCvox.i - vx,
				BCCvox.j - vy,
				BCCvox.k - vz);

			// P1 is the starting BCC point
			P1[0] = vx + vy - vz;
			P1[1] = vx - vy + vz;
			P1[2] = -vx + vy + vz;

			// Sort 
			int alpha_GE_beta = (ga.i >= ga.j);
			int beta_GE_gamma = (ga.j >= ga.k);
			int alpha_GE_gamma = (ga.i >= ga.k);

			I maxParameter = std::max(ga.i, std::max(ga.j, ga.k));
			I minParameter = std::min(ga.i, std::min(ga.j, ga.k));
			I midParameter = ga.i + ga.j + ga.k - maxParameter - minParameter;

			int i = (alpha_GE_beta * 4 + beta_GE_gamma * 2 + alpha_GE_gamma);
			int bm = i, bb = i;
			// We have six cases, the above sort decides which we 
			// use, later on we permute according to this
			I1[0] = (i == 7); I2[0] = (i == 3);
			I1[1] = (i == 5); I2[1] = (i == 2);
			I1[2] = (i == 4); I2[2] = (i == 0);

			// The first two points come for free
			// diagonal
			P2[0] = 1 + P1[0]; 
			P2[1] = 1 + P1[1]; 
			P2[2] = 1 + P1[2]; 
		
			P3[0] = P1[0] + 1 - (I1[2] + I2[2])*2; 
			P3[1] = P1[1] + 1					 - (I2[0] + I2[1])*2; 
			P3[2] = P1[2] - 1 + (I1[2] + I2[2])*2 + (I2[0] + I2[1])*2;  

			P4[0] = P1[0] + 2 - (I1[1] + I1[2])*2 - (I2[1] + I2[2])*2; 
			P4[1] = P1[1] + 0 + (I1[1] + I1[2])*2;  
			P4[2] = P1[2] + 0					 + (I2[1] + I2[2])*2; 


			i = P1[0];
			int j = P1[1];
			int k = P1[2];

			sisl::vector3<I> offset = getSitePosition(i,j,k) - pt;
			I d = offset.norm();

			offset = getSitePosition(P2[0], P2[1], P2[2]) - pt;
			if(offset.norm() < d) {
				d = offset.norm();
				i = P2[0];
				j = P2[1];
				k = P2[2];
			}

			offset = getSitePosition(P3[0], P3[1], P3[2]) - pt;
			if(offset.norm() < d) {
				d = offset.norm();
				i = P3[0];
				j = P3[1];
				k = P3[2];
			}

			offset = getSitePosition(P4[0], P4[1], P4[2]) - pt;
			if(offset.norm() < d) {
				d = offset.norm();
				i = P4[0];
				j = P4[1];
				k = P4[2];
			}
		return vector3<int>(i,j,k);
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
		for(int i = 0; i < res*2; i+=2)
			for(int j = 0; j < res*2; j+=2)
				for(int k = 0; k < res*2; k++){
					int ii = i + (k%2);
					int jj = j + (k%2);

					SV(ii,jj,k, lambda(ii,jj,k));
				}
	}
	
	//
	virtual void spatialFilter(std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> lambda) {
		bcc_odd <GF, O, I> src = *this; // Copy the lattice

		std::function<O(const int &, const int &, const int &)> accessFunc = 
			[&](const int &x, const int &y, const int &z) ->  O {
				return (O) src.GV(x,y,z); 
			};
		
		#pragma omp parallel for
		for(int i = 0; i < res*2; i+=2)
			for(int j = 0; j < res*2; j+=2)
				for(int k = 0; k < res*2; k++){
					int ii = i + (k%2);
					int jj = j + (k%2);
					SV(ii,jj,k, lambda(ii,jj,k, accessFunc));
				}
	} 
	//
	virtual void frequencyFilter(std::function<O(const I &u, const I &v, const I &w, const I &h)> lambda) {
		array3<O> fa(res-1, res-1, res-1);
		array3<O> fb(res, res, res);

		O invdh2 = 1./(dh*dh);

		#pragma omp parallel for
		for(int i = 0; i < res * 2; i+=2)
			for(int j = 0; j < res * 2; j+=2)
				for(int k = 0; k < res; k++){
					if(i != 0 && j != 0 && k != 0)
						fa(i/2 - 1, j/2 - 1, k - 1) = GV(i, j, k*2);
					fb(i/2, j/2, k) = GV(i+1, j+1, k*2 + 1);
				}

		fftw_plan
		p = fftw_plan_r2r_3d(res - 1, res - 1, res - 1,
							fa.getArray(), fa.getArray(), 
							FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, 
							FFTW_ESTIMATE); // TYPE I
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_3d(res, res, res,
							fb.getArray(), fb.getArray(), 
							FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, 
							FFTW_ESTIMATE); // TYPE II
		fftw_execute(p);
		fftw_destroy_plan(p);

		// Construct the transform 'F~' from the above transforms
		#pragma omp parallel for
		for(int i = 0; i < res; i++) 
			for(int j = 0; j < res; j++)
				for(int k = 0; k < res; k++) {
					int _k = res - k - 2;
					O f0 = 0, f0swap = 0;
					if(i < res - 1 && j < res - 1 && k < res - 1)
						f0 = fa(i, j, k);
					L(i, j, k) = (f0 + fb(i, j, k)) * 0.125;
					if(k == res - 1)
						L(i , j, k + res) = 0;
					else {
						if( i < res-1 && j < res - 1 && _k < res - 1)
							f0swap = fa(i, j, _k);
						L(i, j, k + res) = (-f0swap + fb(i,j,_k)) * 0.125; 
					}
				}

		#pragma omp parallel for
		for(int i = 0; i < res; i++)
			for(int j = 0; j < res; j++)
				for(int k = 0; k < res; k++) {
					I u = 2. * M_PI *(1./(4.*res)) * (i+1);
					I v = 2. * M_PI *(1./(4.*res)) * (j+1);
					I wl = 2. * M_PI *(1./(4.*res)) * (k+1);
					I wu = 2. * M_PI *(1./(4.*res)) * (k+1+res);

					L(i,j,k) *= lambda(u,v,wl,dh);
					L(i,j,k+res) *= lambda(u,v,wu,dh);
				}

		#pragma omp parallel for
		for(int i = 0; i < res; i++)
			for(int j = 0; j < res; j++)
				for(int k = 0; k < res; k++) {
						int _k = (res - k - 2);
						if(i < res - 1 && j < res - 1 && k < res - 1) {
							if(_k < 0) 
								fa(i, j, k) = L(i, j, k) - L(i , j, res + res- 1);
							else
								fa(i, j, k) = L(i, j, k) - L(i,  j, _k + res);
						}
						if(_k < 0)
							fb(i,j,k) = L(i,j,k) + L(i, j, res - 1);
						else
							fb(i,j,k) = L(i,j,k) + L(i, j, _k + res);
					}

		// Setup the arrays for the inverse transform
		p = fftw_plan_r2r_3d(res - 1, res - 1, res - 1,
							fa.getArray(), fa.getArray(), 
							FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, 
							FFTW_ESTIMATE); // TYPE I

		fftw_execute(p);
		fftw_destroy_plan(p);

		p = fftw_plan_r2r_3d(res, res, res,
							fb.getArray(), fb.getArray(), 
							FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, 
							FFTW_ESTIMATE); // TYPE II

		fftw_execute(p);
		fftw_destroy_plan(p);

		O scale = 4./O(res*res*res) * 0.125;
		#pragma omp parallel for
		for(int i = 0; i < res*2; i+=2)
			for(int j = 0; j < res*2; j+=2)
				for(int k = 0; k < res; k++) {
					if(i != 0 && j != 0 && k != 0)
						SV(i, j, k*2,
							fa(i/2 - 1, j/2 - 1, k - 1) * scale);
					else
						SV(i, j, k*2, 0);
					SV(i + 1, j + 1, k*2 + 1,
						fb(i/2,j/2,k) * scale);
				}
	} 


	// Convolution of two functions
	//inline cartesian_odd<GF, O, I> operator*(const cartesian_odd<GF, O, I> &v) {throw;};

	// Copy a function
	inline bcc_odd<GF, O, I> operator=(const bcc_odd<GF, O, I> &v) {throw;};

	virtual std::string getLatticeName() const {
		return std::string("bcc");
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

#endif