#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _LINEAR_BCC_BOX_H_
#define _LINEAR_BCC_BOX_H_

namespace sisl{
template <class O, class I>
class linear_bcc_box : public basis_function <O,I> {
public:
	static std::string getBasisName(){
		return std::string("linear");
	}
	static const O M(const I &xi, const I &yi, const I &zi) { 
		O x = std::max(fabs(xi),std::max(fabs(yi),fabs(zi)));
		O z = std::min(fabs(xi),std::min(fabs(yi),fabs(zi)));
		O y = fabs(xi)+fabs(yi)+fabs(zi)-x-z;
		
		if(x+y > 0) return 0;
		return (2.-x-y)*(1./8.);
	}
	static const O M(const vector3<I> &p) { return linear_bcc_box<O,I>::M(p.i, p.j, p.k); }

	// This function should return the intersection of the closure of the support of
	// the generator and the lattice.
	static std::vector<std::tuple<int,int,int>> getSupport() {
		using namespace std;
		return {make_tuple(-2, -2, -2), make_tuple(-2, -2, 0), make_tuple(-2, -2, 2), make_tuple(-2, 0, -2), make_tuple(-2, 0, 0), make_tuple(-2, 0, 2), make_tuple(-2, 2, -2), make_tuple(-2, 2, 0), make_tuple(-2, 2, 2), make_tuple(-1, -1, -1), make_tuple(-1, -1, 1), make_tuple(-1, 1, -1), make_tuple(-1, 1, 1), make_tuple(0, -2, -2), make_tuple(0, -2, 0), make_tuple(0, -2, 2), make_tuple(0, 0, -2), make_tuple(0, 0, 0), make_tuple(0, 0, 2), make_tuple(0, 2, -2), make_tuple(0, 2, 0), make_tuple(0, 2, 2), make_tuple(1, -1, -1), make_tuple(1, -1, 1), make_tuple(1, 1, -1), make_tuple(1, 1, 1), make_tuple(2, -2, -2), make_tuple(2, -2, 0), make_tuple(2, -2, 2), make_tuple(2, 0, -2), make_tuple(2, 0, 0), make_tuple(2, 0, 2), make_tuple(2, 2, -2), make_tuple(2, 2, 0), make_tuple(2, 2, 2)};
	}
	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int>> getEffectiveSupport(const vector3<I> &p) {
		using namespace std;
		return {make_tuple(-2, -2, -2), make_tuple(-2, -2, 0), make_tuple(-2, -2, 2), make_tuple(-2, 0, -2), make_tuple(-2, 0, 0), make_tuple(-2, 0, 2), make_tuple(-2, 2, -2), make_tuple(-2, 2, 0), make_tuple(-2, 2, 2), make_tuple(-1, -1, -1), make_tuple(-1, -1, 1), make_tuple(-1, 1, -1), make_tuple(-1, 1, 1), make_tuple(0, -2, -2), make_tuple(0, -2, 0), make_tuple(0, -2, 2), make_tuple(0, 0, -2), make_tuple(0, 0, 0), make_tuple(0, 0, 2), make_tuple(0, 2, -2), make_tuple(0, 2, 0), make_tuple(0, 2, 2), make_tuple(1, -1, -1), make_tuple(1, -1, 1), make_tuple(1, 1, -1), make_tuple(1, 1, 1), make_tuple(2, -2, -2), make_tuple(2, -2, 0), make_tuple(2, -2, 2), make_tuple(2, 0, -2), make_tuple(2, 0, 0), make_tuple(2, 0, 2), make_tuple(2, 2, -2), make_tuple(2, 2, 0), make_tuple(2, 2, 2)};
	}
	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi2Norm(){throw "basis_function()::getBeppoLevi2Norm() - Does not exist for linear box spline";};
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> autoCorrelation(){throw "basis_function()::autoCorrelation() - Not Implemented!";};

	// This does the actual semi-descrete convolution sum
	// the idea behind this, is to allow basis functions
	// to provide a potentially optimized version of
	// the convolution sum for particular lattices.
	static const O convolutionSum(const vector3<I> &p, const shift_invariant_space3<linear_bcc_box, O, I> *L) {

		int P1[3],P2[3],P3[3],P4[3];
		int Q1[3],Q2[3],Q3[3],Q4[3];
		int I1[3],I2[3];

		O r = 0;
		I h = L->getScale();

		I x = p.i/h;
		I y = p.j/h;
		I z = p.k/h;

		vector3<I> BCCvox(
			(x + y) / 2,
			(x + z) / 2,
			(y + z) / 2);

		int vx = (int)floor(BCCvox.i),
			vy = (int)floor(BCCvox.j),
			vz = (int)floor(BCCvox.k);

		vector3<O> ga(
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

		// 
		I1[0] = (i == 7); I2[0] = (i == 3);
		I1[1] = (i == 5); I2[1] = (i == 2);
		I1[2] = (i == 4); I2[2] = (i == 0);

		// The first two points come for free
		// diagonal
		P2[0] = 1 + P1[0];
		P2[1] = 1 + P1[1];
		P2[2] = 1 + P1[2]; 
	
		// The rest are determinied by the sort
		P3[0] = P1[0] + 1 - (I1[2] + I2[2])*2; 
		P3[1] = P1[1] + 1                     - (I2[0] + I2[1])*2; 
		P3[2] = P1[2] - 1 + (I1[2] + I2[2])*2 + (I2[0] + I2[1])*2;  

		P4[0] = P1[0] + 2 - (I1[1] + I1[2])*2 - (I2[1] + I2[2])*2; 
		P4[1] = P1[1] + 0 + (I1[1] + I1[2])*2;
		P4[2] = P1[2] + 0                     + (I2[1] + I2[2])*2;

		O D1 = L->GV(P1[0], P1[1], P1[2]) * (1-maxParameter);
		O D2 = L->GV(P2[0], P2[1], P2[2]) * minParameter;
		O D3 = L->GV(P3[0], P3[1], P3[2]) * (maxParameter-midParameter);
		O D4 = L->GV(P4[0], P4[1], P4[2]) * (midParameter-minParameter);
		return D1 + D2 + D3 + D4;
	}

	//Filters	
	static O autocorrelationFilter(const I &u, const I &v, const I &w, const I &h) {
		return 1./((1./15.) * (6. + 6.*cos(u)*cos(v)*cos(w) + cos(2.*u) + cos(2.*v) + cos(2.*w)));
	}

	static O poissonFilter(const I &u, const I &v, const I &w, const I &h) {
		O invh = 1./(h*h);
		return 1./(invh*(1./6.) * (-15. + 16.*cos(u)*cos(v)*cos(w) - cos(2.*u)*cos(2.*v)*cos(2.*w)));
	}
};

};
#endif