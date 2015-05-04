#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _QUNTIC_BOX_H_
#define _QUNTIC_BOX_H_

namespace sisl{
template <class O, class I>
class quintic_box : public basis_function <O,I> {
public:
	static std::string getBasisName(){
		return std::string("quintic");
	}
	static const O M(const I &xi, const I &yi, const I &zi) { 
		O x = std::max(fabs(xi),std::max(fabs(yi),fabs(zi)));
		O z = std::min(fabs(xi),std::min(fabs(yi),fabs(zi)));
		O y = fabs(xi)+fabs(yi)+fabs(zi)-x-z;
		
		if((x+y) > 4) return 0;

		O xy4 = (x + y - 4.);
		O xz2 = (x + z - 2.);
		O yz2 = (y + z - 2.);
		O xy2 = (x + y - 2.);

		O ret = 0;
		xy4 = xy4*xy4*xy4;
		xz2 = xz2*xz2*xz2;
		yz2 = yz2*yz2*yz2;
		xy2 = xy2*xy2*xy2;

		if( (x+y) < 2) { // Region R1
			ret = (1./3840.) * xy4 * (-3.*x*y -5.*z*z +2.*x +2.*y +20.*z +x*x +y*y -24.);
			ret += (1./1920.) * xz2 * (x*x -9.*x -3.*x*z +10.*y -5.*y*y +14. +11.*z +z*z);
			ret += (1./1920.) * yz2 * (46. -30.*x -z-y +3.*z*y +5.*x*x -y*y -z*z);
			ret -= (1./960.)* xy2 * (x*x +x -3.*x*y -5.*z*z +y*y +y -6.);
			return ret * 4;
		}else if(x+z < 2){ //Region R2
			ret =  (1./3840.) * xy4 * (-3.*x*y - 5.*z*z + 2.*x + 2.*y + 20.*z + x*x + y*y - 24.) -
					(1./1920.)  * xz2 * (-z*z - 11.*z + 3.*x*z - 14. + 5.*y*y + 9.*x - 10.*y - x*x) -
					(1./1920.)  * yz2 * (-46. + z + 30.*x + y - 3.*z*y - 5.*x*x + y*y + z*z);
			return ret * 4;
		}else if(y+z < 2){ // Region R3
			if((x-z) > 2){ // Region R3,A
				ret = (1./3840.) * xy4 * ( -x*x +8.*x +3.*x*y -y*y +5.*z*z -16. -12.*y);
				return ret * 4;
			}else{ // Region R3,B
				ret = xy4 * (-3.*x*y -5.*z*z +2.*x +2.*y +20.*z +x*x +y*y -24.)/3840. -
						yz2 * (30.*x +z -46. -3.*y*z +y -5.*x*x +y*y+z*z)/1920.;
				return ret * 4;
			}
		}
		// R4
		ret = xy4 * (-3.*x*y -5.*z*z +2.*x +2.*y +20.*z +x*x +y*y -24.)/3840.;
		return ret * 4;
	}
	static const O M(const vector3<I> &p) { return quintic_box<O,I>::M(p.i, p.j, p.k); }

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
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi2Norm(){
		using namespace std;
		return {make_tuple(-6, 0, 0, 1./20160.), make_tuple(-5, -1, -1, 23./24192.), make_tuple(-5, -1, 1, 23./24192.), make_tuple(-5, 1, -1, 23./24192.), make_tuple(-5, 1, 1, 23./24192.), make_tuple(-4, -2, -2, 61./30240.), make_tuple(-4, -2, 0, 361./60480.), make_tuple(-4, -2, 2, 61./30240.), make_tuple(-4, 0, -2, 361./60480.), make_tuple(-4, 0, 0, 13./1512.), make_tuple(-4, 0, 2, 361./60480.), make_tuple(-4, 2, -2, 61./30240.), make_tuple(-4, 2, 0, 361./60480.), make_tuple(-4, 2, 2, 61./30240.), make_tuple(-3, -3, -3, 11./24192.), make_tuple(-3, -3, -1, 1019./120960.), make_tuple(-3, -3, 1, 1019./120960.), make_tuple(-3, -3, 3, 11./24192.), make_tuple(-3, -1, -3, 1019./120960.), make_tuple(-3, -1, -1, -(487./60480.)), make_tuple(-3, -1, 1, -(487./60480.)), make_tuple(-3, -1, 3, 1019./120960.), make_tuple(-3, 1, -3, 1019./120960.), make_tuple(-3, 1, -1, -(487./60480.)), make_tuple(-3, 1, 1, -(487./60480.)), make_tuple(-3, 1, 3, 1019./120960.), make_tuple(-3, 3, -3, 11./24192.), make_tuple(-3, 3, -1, 1019./120960.), make_tuple(-3, 3, 1, 1019./120960.), make_tuple(-3, 3, 3, 11./24192.), make_tuple(-2, -4, -2, 61./30240.), make_tuple(-2, -4, 0, 361./60480.), make_tuple(-2, -4, 2, 61./30240.), make_tuple(-2, -2, -4, 61./30240.), make_tuple(-2, -2, -2, 1./252.), make_tuple(-2, -2, 0, -(7./120.)), make_tuple(-2, -2, 2, 1./252.), make_tuple(-2, -2, 4, 61./30240.), make_tuple(-2, 0, -4, 361./60480.), make_tuple(-2, 0, -2, -(7./120.)), make_tuple(-2, 0, 0, -(3071./60480.)), make_tuple(-2, 0, 2, -(7./120.)), make_tuple(-2, 0, 4, 361./60480.), make_tuple(-2, 2, -4, 61./30240.), make_tuple(-2, 2, -2, 1./252.), make_tuple(-2, 2, 0, -(7./120.)), make_tuple(-2, 2, 2, 1./252.), make_tuple(-2, 2, 4, 61./30240.), make_tuple(-2, 4, -2, 61./30240.), make_tuple(-2, 4, 0, 361./60480.), make_tuple(-2, 4, 2, 61./30240.), make_tuple(-1, -5, -1, 23./24192.), make_tuple(-1, -5, 1, 23./24192.), make_tuple(-1, -3, -3, 1019./120960.), make_tuple(-1, -3, -1, -(487./60480.)), make_tuple(-1, -3, 1, -(487./60480.)), make_tuple(-1, -3, 3, 1019./120960.), make_tuple(-1, -1, -5, 23./24192.), make_tuple(-1, -1, -3, -(487./60480.)), make_tuple(-1, -1, -1, -(479./13440.)), make_tuple(-1, -1, 1, -(479./13440.)), make_tuple(-1, -1, 3, -(487./60480.)), make_tuple(-1, -1, 5, 23./24192.), make_tuple(-1, 1, -5, 23./24192.), make_tuple(-1, 1, -3, -(487./60480.)), make_tuple(-1, 1, -1, -(479./13440.)), make_tuple(-1, 1, 1, -(479./13440.)), make_tuple(-1, 1, 3, -(487./60480.)), make_tuple(-1, 1, 5, 23./24192.), make_tuple(-1, 3, -3, 1019./120960.), make_tuple(-1, 3, -1, -(487./60480.)), make_tuple(-1, 3, 1, -(487./60480.)), make_tuple(-1, 3, 3, 1019./120960.), make_tuple(-1, 5, -1, 23./24192.), make_tuple(-1, 5, 1, 23./24192.), make_tuple(0, -6, 0, 1./20160.), make_tuple(0, -4, -2, 361./60480.), make_tuple(0, -4, 0, 13./1512.), make_tuple(0, -4, 2, 361./60480.), make_tuple(0, -2, -4, 361./60480.), make_tuple(0, -2, -2, -(7./120.)), make_tuple(0, -2, 0, -(3071./60480.)), make_tuple(0, -2, 2, -(7./120.)), make_tuple(0, -2, 4, 361./60480.), make_tuple(0, 0, -6, 1./20160.), make_tuple(0, 0, -4, 13./1512.), make_tuple(0, 0, -2, -(3071./60480.)), make_tuple(0, 0, 0, 3701./3780.), make_tuple(0, 0, 2, -(3071./60480.)), make_tuple(0, 0, 4, 13./1512.), make_tuple(0, 0, 6, 1./20160.), make_tuple(0, 2, -4, 361./60480.), make_tuple(0, 2, -2, -(7./120.)), make_tuple(0, 2, 0, -(3071./60480.)), make_tuple(0, 2, 2, -(7./120.)), make_tuple(0, 2, 4, 361./60480.), make_tuple(0, 4, -2, 361./60480.), make_tuple(0, 4, 0, 13./1512.), make_tuple(0, 4, 2, 361./60480.), make_tuple(0, 6, 0, 1./20160.), make_tuple(1, -5, -1, 23./24192.), make_tuple(1, -5, 1, 23./24192.), make_tuple(1, -3, -3, 1019./120960.), make_tuple(1, -3, -1, -(487./60480.)), make_tuple(1, -3, 1, -(487./60480.)), make_tuple(1, -3, 3, 1019./120960.), make_tuple(1, -1, -5, 23./24192.), make_tuple(1, -1, -3, -(487./60480.)), make_tuple(1, -1, -1, -(479./13440.)), make_tuple(1, -1, 1, -(479./13440.)), make_tuple(1, -1, 3, -(487./60480.)), make_tuple(1, -1, 5, 23./24192.), make_tuple(1, 1, -5, 23./24192.), make_tuple(1, 1, -3, -(487./60480.)), make_tuple(1, 1, -1, -(479./13440.)), make_tuple(1, 1, 1, -(479./13440.)), make_tuple(1, 1, 3, -(487./60480.)), make_tuple(1, 1, 5, 23./24192.), make_tuple(1, 3, -3, 1019./120960.), make_tuple(1, 3, -1, -(487./60480.)), make_tuple(1, 3, 1, -(487./60480.)), make_tuple(1, 3, 3, 1019./120960.), make_tuple(1, 5, -1, 23./24192.), make_tuple(1, 5, 1, 23./24192.), make_tuple(2, -4, -2, 61./30240.), make_tuple(2, -4, 0, 361./60480.), make_tuple(2, -4, 2, 61./30240.), make_tuple(2, -2, -4, 61./30240.), make_tuple(2, -2, -2, 1./252.), make_tuple(2, -2, 0, -(7./120.)), make_tuple(2, -2, 2, 1./252.), make_tuple(2, -2, 4, 61./30240.), make_tuple(2, 0, -4, 361./60480.), make_tuple(2, 0, -2, -(7./120.)), make_tuple(2, 0, 0, -(3071./60480.)), make_tuple(2, 0, 2, -(7./120.)), make_tuple(2, 0, 4, 361./60480.), make_tuple(2, 2, -4, 61./30240.), make_tuple(2, 2, -2, 1./252.), make_tuple(2, 2, 0, -(7./120.)), make_tuple(2, 2, 2, 1./252.), make_tuple(2, 2, 4, 61./30240.), make_tuple(2, 4, -2, 61./30240.), make_tuple(2, 4, 0, 361./60480.), make_tuple(2, 4, 2, 61./30240.), make_tuple(3, -3, -3, 11./24192.), make_tuple(3, -3, -1, 1019./120960.), make_tuple(3, -3, 1, 1019./120960.), make_tuple(3, -3, 3, 11./24192.), make_tuple(3, -1, -3, 1019./120960.), make_tuple(3, -1, -1, -(487./60480.)), make_tuple(3, -1, 1, -(487./60480.)), make_tuple(3, -1, 3, 1019./120960.), make_tuple(3, 1, -3, 1019./120960.), make_tuple(3, 1, -1, -(487./60480.)), make_tuple(3, 1, 1, -(487./60480.)), make_tuple(3, 1, 3, 1019./120960.), make_tuple(3, 3, -3, 11./24192.), make_tuple(3, 3, -1, 1019./120960.), make_tuple(3, 3, 1, 1019./120960.), make_tuple(3, 3, 3, 11./24192.), make_tuple(4, -2, -2, 61./30240.), make_tuple(4, -2, 0, 361./60480.), make_tuple(4, -2, 2, 61./30240.), make_tuple(4, 0, -2, 361./60480.), make_tuple(4, 0, 0, 13./1512.), make_tuple(4, 0, 2, 361./60480.), make_tuple(4, 2, -2, 61./30240.), make_tuple(4, 2, 0, 361./60480.), make_tuple(4, 2, 2, 61./30240.), make_tuple(5, -1, -1, 23./24192.), make_tuple(5, -1, 1, 23./24192.), make_tuple(5, 1, -1, 23./24192.), make_tuple(5, 1, 1, 23./24192.), make_tuple(6, 0, 0, 1./20160.)};
	};
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> autoCorrelation(){
		using namespace std;
		return {make_tuple(-6, 0, 0, 6.4236175347286453e-08), make_tuple(-5, -1, -1, 2.1454882565993675e-06), make_tuple(-5, -1, 1, 2.1454882565993675e-06), make_tuple(-5, 1, -1, 2.1454882565993675e-06), make_tuple(-5, 1, 1, 2.1454882565993675e-06), make_tuple(-4, -2, -2, 3.0319474763919209e-06), make_tuple(-4, -2, 0, 4.2331639553861779e-05), make_tuple(-4, -2, 2, 3.0319474763919209e-06), make_tuple(-4, 0, -2, 4.2331639553861779e-05), make_tuple(-4, 0, 0, 4.8983937872826766e-04), make_tuple(-4, 0, 2, 4.2331639553861779e-05), make_tuple(-4, 2, -2, 3.0319474763919209e-06), make_tuple(-4, 2, 0, 4.2331639553861779e-05), make_tuple(-4, 2, 2, 3.0319474763919209e-06), make_tuple(-3, -3, -3, 2.6979193645860314e-07), make_tuple(-3, -3, -1, 4.6185810074698963e-05), make_tuple(-3, -3, 1, 4.6185810074698963e-05), make_tuple(-3, -3, 3, 2.6979193645860314e-07), make_tuple(-3, -1, -3, 4.6185810074698963e-05), make_tuple(-3, -1, -1, 2.7537534481978927e-03), make_tuple(-3, -1, 1, 2.7537534481978927e-03), make_tuple(-3, -1, 3, 4.6185810074698963e-05), make_tuple(-3, 1, -3, 4.6185810074698963e-05), make_tuple(-3, 1, -1, 2.7537534481978927e-03), make_tuple(-3, 1, 1, 2.7537534481978927e-03), make_tuple(-3, 1, 3, 4.6185810074698963e-05), make_tuple(-3, 3, -3, 2.6979193645860314e-07), make_tuple(-3, 3, -1, 4.6185810074698963e-05), make_tuple(-3, 3, 1, 4.6185810074698963e-05), make_tuple(-3, 3, 3, 2.6979193645860314e-07), make_tuple(-2, -4, -2, 3.0319474763919209e-06), make_tuple(-2, -4, 0, 4.2331639553861779e-05), make_tuple(-2, -4, 2, 3.0319474763919209e-06), make_tuple(-2, -2, -4, 3.0319474763919209e-06), make_tuple(-2, -2, -2, 1.4920264920264921e-03), make_tuple(-2, -2, 0, 8.4620626287292954e-03), make_tuple(-2, -2, 2, 1.4920264920264921e-03), make_tuple(-2, -2, 4, 3.0319474763919209e-06), make_tuple(-2, 0, -4, 4.2331639553861779e-05), make_tuple(-2, 0, -2, 8.4620626287292954e-03), make_tuple(-2, 0, 0, 3.8925978856534413e-02), make_tuple(-2, 0, 2, 8.4620626287292954e-03), make_tuple(-2, 0, 4, 4.2331639553861779e-05), make_tuple(-2, 2, -4, 3.0319474763919209e-06), make_tuple(-2, 2, -2, 1.4920264920264921e-03), make_tuple(-2, 2, 0, 8.4620626287292954e-03), make_tuple(-2, 2, 2, 1.4920264920264921e-03), make_tuple(-2, 2, 4, 3.0319474763919209e-06), make_tuple(-2, 4, -2, 3.0319474763919209e-06), make_tuple(-2, 4, 0, 4.2331639553861779e-05), make_tuple(-2, 4, 2, 3.0319474763919209e-06), make_tuple(-1, -5, -1, 2.1454882565993675e-06), make_tuple(-1, -5, 1, 2.1454882565993675e-06), make_tuple(-1, -3, -3, 4.6185810074698963e-05), make_tuple(-1, -3, -1, 2.7537534481978927e-03), make_tuple(-1, -3, 1, 2.7537534481978927e-03), make_tuple(-1, -3, 3, 4.6185810074698963e-05), make_tuple(-1, -1, -5, 2.1454882565993675e-06), make_tuple(-1, -1, -3, 2.7537534481978927e-03), make_tuple(-1, -1, -1, 5.3815821524154858e-02), make_tuple(-1, -1, 1, 5.3815821524154858e-02), make_tuple(-1, -1, 3, 2.7537534481978927e-03), make_tuple(-1, -1, 5, 2.1454882565993675e-06), make_tuple(-1, 1, -5, 2.1454882565993675e-06), make_tuple(-1, 1, -3, 2.7537534481978927e-03), make_tuple(-1, 1, -1, 5.3815821524154858e-02), make_tuple(-1, 1, 1, 5.3815821524154858e-02), make_tuple(-1, 1, 3, 2.7537534481978927e-03), make_tuple(-1, 1, 5, 2.1454882565993675e-06), make_tuple(-1, 3, -3, 4.6185810074698963e-05), make_tuple(-1, 3, -1, 2.7537534481978927e-03), make_tuple(-1, 3, 1, 2.7537534481978927e-03), make_tuple(-1, 3, 3, 4.6185810074698963e-05), make_tuple(-1, 5, -1, 2.1454882565993675e-06), make_tuple(-1, 5, 1, 2.1454882565993675e-06), make_tuple(0, -6, 0, 6.4236175347286453e-08), make_tuple(0, -4, -2, 4.2331639553861779e-05), make_tuple(0, -4, 0, 4.8983937872826766e-04), make_tuple(0, -4, 2, 4.2331639553861779e-05), make_tuple(0, -2, -4, 4.2331639553861779e-05), make_tuple(0, -2, -2, 8.4620626287292954e-03), make_tuple(0, -2, 0, 3.8925978856534413e-02), make_tuple(0, -2, 2, 8.4620626287292954e-03), make_tuple(0, -2, 4, 4.2331639553861779e-05), make_tuple(0, 0, -6, 6.4236175347286453e-08), make_tuple(0, 0, -4, 4.8983937872826766e-04), make_tuple(0, 0, -2, 3.8925978856534413e-02), make_tuple(0, 0, 0, 1.5115625115625114e-01), make_tuple(0, 0, 2, 3.8925978856534413e-02), make_tuple(0, 0, 4, 4.8983937872826766e-04), make_tuple(0, 0, 6, 6.4236175347286453e-08), make_tuple(0, 2, -4, 4.2331639553861779e-05), make_tuple(0, 2, -2, 8.4620626287292954e-03), make_tuple(0, 2, 0, 3.8925978856534413e-02), make_tuple(0, 2, 2, 8.4620626287292954e-03), make_tuple(0, 2, 4, 4.2331639553861779e-05), make_tuple(0, 4, -2, 4.2331639553861779e-05), make_tuple(0, 4, 0, 4.8983937872826766e-04), make_tuple(0, 4, 2, 4.2331639553861779e-05), make_tuple(0, 6, 0, 6.4236175347286453e-08), make_tuple(1, -5, -1, 2.1454882565993675e-06), make_tuple(1, -5, 1, 2.1454882565993675e-06), make_tuple(1, -3, -3, 4.6185810074698963e-05), make_tuple(1, -3, -1, 2.7537534481978927e-03), make_tuple(1, -3, 1, 2.7537534481978927e-03), make_tuple(1, -3, 3, 4.6185810074698963e-05), make_tuple(1, -1, -5, 2.1454882565993675e-06), make_tuple(1, -1, -3, 2.7537534481978927e-03), make_tuple(1, -1, -1, 5.3815821524154858e-02), make_tuple(1, -1, 1, 5.3815821524154858e-02), make_tuple(1, -1, 3, 2.7537534481978927e-03), make_tuple(1, -1, 5, 2.1454882565993675e-06), make_tuple(1, 1, -5, 2.1454882565993675e-06), make_tuple(1, 1, -3, 2.7537534481978927e-03), make_tuple(1, 1, -1, 5.3815821524154858e-02), make_tuple(1, 1, 1, 5.3815821524154858e-02), make_tuple(1, 1, 3, 2.7537534481978927e-03), make_tuple(1, 1, 5, 2.1454882565993675e-06), make_tuple(1, 3, -3, 4.6185810074698963e-05), make_tuple(1, 3, -1, 2.7537534481978927e-03), make_tuple(1, 3, 1, 2.7537534481978927e-03), make_tuple(1, 3, 3, 4.6185810074698963e-05), make_tuple(1, 5, -1, 2.1454882565993675e-06), make_tuple(1, 5, 1, 2.1454882565993675e-06), make_tuple(2, -4, -2, 3.0319474763919209e-06), make_tuple(2, -4, 0, 4.2331639553861779e-05), make_tuple(2, -4, 2, 3.0319474763919209e-06), make_tuple(2, -2, -4, 3.0319474763919209e-06), make_tuple(2, -2, -2, 1.4920264920264921e-03), make_tuple(2, -2, 0, 8.4620626287292954e-03), make_tuple(2, -2, 2, 1.4920264920264921e-03), make_tuple(2, -2, 4, 3.0319474763919209e-06), make_tuple(2, 0, -4, 4.2331639553861779e-05), make_tuple(2, 0, -2, 8.4620626287292954e-03), make_tuple(2, 0, 0, 3.8925978856534413e-02), make_tuple(2, 0, 2, 8.4620626287292954e-03), make_tuple(2, 0, 4, 4.2331639553861779e-05), make_tuple(2, 2, -4, 3.0319474763919209e-06), make_tuple(2, 2, -2, 1.4920264920264921e-03), make_tuple(2, 2, 0, 8.4620626287292954e-03), make_tuple(2, 2, 2, 1.4920264920264921e-03), make_tuple(2, 2, 4, 3.0319474763919209e-06), make_tuple(2, 4, -2, 3.0319474763919209e-06), make_tuple(2, 4, 0, 4.2331639553861779e-05), make_tuple(2, 4, 2, 3.0319474763919209e-06), make_tuple(3, -3, -3, 2.6979193645860314e-07), make_tuple(3, -3, -1, 4.6185810074698963e-05), make_tuple(3, -3, 1, 4.6185810074698963e-05), make_tuple(3, -3, 3, 2.6979193645860314e-07), make_tuple(3, -1, -3, 4.6185810074698963e-05), make_tuple(3, -1, -1, 2.7537534481978927e-03), make_tuple(3, -1, 1, 2.7537534481978927e-03), make_tuple(3, -1, 3, 4.6185810074698963e-05), make_tuple(3, 1, -3, 4.6185810074698963e-05), make_tuple(3, 1, -1, 2.7537534481978927e-03), make_tuple(3, 1, 1, 2.7537534481978927e-03), make_tuple(3, 1, 3, 4.6185810074698963e-05), make_tuple(3, 3, -3, 2.6979193645860314e-07), make_tuple(3, 3, -1, 4.6185810074698963e-05), make_tuple(3, 3, 1, 4.6185810074698963e-05), make_tuple(3, 3, 3, 2.6979193645860314e-07), make_tuple(4, -2, -2, 3.0319474763919209e-06), make_tuple(4, -2, 0, 4.2331639553861779e-05), make_tuple(4, -2, 2, 3.0319474763919209e-06), make_tuple(4, 0, -2, 4.2331639553861779e-05), make_tuple(4, 0, 0, 4.8983937872826766e-04), make_tuple(4, 0, 2, 4.2331639553861779e-05), make_tuple(4, 2, -2, 3.0319474763919209e-06), make_tuple(4, 2, 0, 4.2331639553861779e-05), make_tuple(4, 2, 2, 3.0319474763919209e-06), make_tuple(5, -1, -1, 2.1454882565993675e-06), make_tuple(5, -1, 1, 2.1454882565993675e-06), make_tuple(5, 1, -1, 2.1454882565993675e-06), make_tuple(5, 1, 1, 2.1454882565993675e-06), make_tuple(6, 0, 0, 6.4236175347286453e-08)};
	};

	inline static const vector3<O> approximateGradient(const int &i, const int &j, const int &k, const shift_invariant_space3< quintic_box, O, I> *L) {
		O sx = O(1./(L->getScale()));
#define W1 (1./6.)
#define W2 (1./48.)
		O dx = (
			(L->GV(i+(1),j+(1),k+(1))  +  L->GV(i+(1),j+(1),k+(-1))  + L->GV(i+(1),j+(-1),k+(1))  + L->GV(i+(1),j+(-1),k+(-1)) )*(-W1) +				
			(L->GV(i+(-1),j+(1),k+(1))  + L->GV(i+(-1),j+(1),k+(-1)) + L->GV(i+(-1),j+(-1),k+(1)) + L->GV(i+(-1),j+(-1),k+(-1)) )*(W1) +
			(L->GV(i+(2),j+(2),k+(2))  + L->GV(i+(2),j+(2),k+(-2)) + L->GV(i+(2),j+(-2),k+(2)) + L->GV(i+(2),j+(-2),k+(-2)) )*(W2) +
			(L->GV(i+(-2),j+(2),k+(2))  + L->GV(i+(-2),j+(2),k+(-2)) + L->GV(i+(-2),j+(-2),k+(2)) + L->GV(i+(-2),j+(-2),k+(-2)) )*(-W2)) * sx;
		
		O dy = (
			(L->GV(i+(1),j+(1),k+(1))  +  L->GV(i+(1),j+(1),k+(-1))  + L->GV(i+(-1),j+(1),k+(1))  + L->GV(i+(-1),j+(1),k+(-1)) )*(-W1) +				
			(L->GV(i+(1),j+(-1),k+(1))  + L->GV(i+(1),j+(-1),k+(-1)) + L->GV(i+(-1),j+(-1),k+(1)) + L->GV(i+(-1),j+(-1),k+(-1)) )*(W1) +
			(L->GV(i+(2),j+(2),k+(2))  + L->GV(i+(2),j+(2),k+(-2)) + L->GV(i+(-2),j+(2),k+(2)) + L->GV(i+(-2),j+(2),k+(-2)) )*(W2) +
			(L->GV(i+(2),j+(-2),k+(2))  + L->GV(i+(2),j+(-2),k+(-2)) + L->GV(i+(-2),j+(-2),k+(2)) + L->GV(i+(-2),j+(-2),k+(-2)) )*(-W2))* sx;

		O dz = (
			(L->GV(i+(1),j+(1),k+(1))  +  L->GV(i+(-1),j+(1),k+(1))  + L->GV(i+(1),j+(-1),k+(1))  + L->GV(i+(-1),j+(-1),k+(1)) )*(-W1) +				
			(L->GV(i+(1),j+(1),k+(-1))  + L->GV(i+(-1),j+(1),k+(-1)) + L->GV(i+(1),j+(-1),k+(-1)) + L->GV(i+(-1),j+(-1),k+(-1)) )*(W1) +
			(L->GV(i+(2),j+(2),k+(2))  + L->GV(i+(-2),j+(2),k+(2)) + L->GV(i+(2),j+(-2),k+(2)) + L->GV(i+(-2),j+(-2),k+(2)) )*(W2) +
			(L->GV(i+(2),j+(2),k+(-2))  + L->GV(i+(-2),j+(2),k+(-2)) + L->GV(i+(2),j+(-2),k+(-2)) + L->GV(i+(-2),j+(-2),k+(-2)) )*(-W2)) * sx ;
#undef W1
#undef W2
		return vector3<O>(dx,dy,dz);
	}

	// This does the actual semi-descrete convolution sum
	// the idea behind this, is to allow basis functions
	// to provide a potentially optimized version of
	// the convolution sum for particular lattices.
	static const O convolutionSum(const vector3<I> &p, const shift_invariant_space3<quintic_box, O, I> *L) {
//
#define RHO_FAST(a, b, g, a2, a3) \
    -a3*(g*b-0.5*a*(g+b)+0.3*a2);

#define RHO22(A, B, G) ( -((A)*(A))*(A) \
    * ((1./6.)*(G)*(B) - (1./12.) *(A) \
    * ((G)+(B))+0.05 * ((A)*(A))))

#define FILL_PPIPED(p0, p1, p2, p3, p4, p5, p6, p7, Q) { \
	p0 = L->GV(x0,        y0,        z0      ); \
	p1 = L->GV(x0 - Q[0], y0 + Q[1], z0 + Q[2]); \
	p2 = L->GV(x0 + Q[0], y0 - Q[1], z0 + Q[2]); \
	p3 = L->GV(x0 + Q[0], y0 + Q[1], z0 - Q[2]); \
	p4 = L->GV(x0 + 2*Q[0], y0,          z0         ); \
	p5 = L->GV(x0,          y0 + 2*Q[1], z0         ); \
	p6 = L->GV(x0,          y0,          z0 + 2*Q[2]); \
	p7 = L->GV(x0 + Q[0],   y0 + Q[1],   z0 + Q[2]   ); \
	}

#define CONV_PPIPED(value, pa, pb, alpha, beta, gamma) { \
    O p123 = pa[1] + pa[2] + pa[3]; \
    O p0 = pa[0] * 4.0; \
    O alpha2 = alpha * alpha; \
    O alpha3 = alpha2 * alpha * (1./6.); \
    O base_rho = RHO_FAST(alpha, beta, gamma, alpha2, alpha3); \
    O base_rho2 = base_rho - 0.5 * alpha3 * alpha; \
    value += (-2.5*p0+4.0*(p123)-2.0*(pb[0]+pb[1]+pb[2]) + pb[3]) * base_rho; \
    alpha2 = base_rho2 + alpha3 * beta; \
    base_rho2 += alpha3 * gamma; \
    value += (p0-2.0*(p123-pa[1])+pb[0]) * (alpha2) \
    + (p0-2.0*(p123-pa[2])+pb[1]) * (base_rho2) \
    + (-.5*p0+pa[3]) \
    * (base_rho2 + alpha2 - alpha3 - base_rho); \
    alpha -= 1.0; alpha2 = beta * beta; \
    alpha3 = (1./6.) * alpha2 * beta; \
    base_rho = RHO_FAST(beta, gamma, alpha, alpha2, alpha3); \
    value += (p0-2.0*(p123-pa[3])+pb[2]) * (base_rho); \
    p0 = -.5*p0; \
    value += (p0+pa[2]) * (base_rho + alpha3 * (alpha - .5*beta)) \
    + (p0 + pa[1]) * RHO22(gamma, alpha, beta-1.0) + (-0.5 * p0) \
    * RHO22(alpha, beta-1.0, gamma-1.0); \
	}

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

		// We have six cases, the above sort decides which we 
		// use, later on we permute according to this
		I1[0] = (i == 7); I2[0] = (i == 3);
		I1[1] = (i == 5); I2[1] = (i == 2);
		I1[2] = (i == 4); I2[2] = (i == 0);

		// The first two points come for free
		// diagonal
		P2[0] = 1; Q1[0] = -1;
		P2[1] = 1; Q1[1] = -1;
		P2[2] = 1; Q1[2] = -1;
	
		// The rest are determinied by the sort
		Q3[0] = P3[0] =  1 - (I1[2] + I2[2])*2; 
		Q3[1] = P3[1] =  1                     - (I2[0] + I2[1])*2; 
		Q3[2] = P3[2] = -1 + (I1[2] + I2[2])*2 + (I2[0] + I2[1])*2;  

		P4[0] =  2 - (I1[1] + I1[2])*2 - (I2[1] + I2[2])*2; 
		P4[1] =  0 + (I1[1] + I1[2])*2;  
		P4[2] =  0                     + (I2[1] + I2[2])*2; 

		Q2[0] = -1 + (I1[1] + I1[2])*2 + (I2[1] + I2[2])*2; 
		Q2[1] =  1 - (I1[1] + I1[2])*2;  
		Q2[2] =  1                     - (I2[1] + I2[2])*2; 

		Q4[0] =  1 - (I1[1] + I2[1])*2 ; 
		Q4[1] = -1 + (I1[1] + I2[1])*2 + (I1[2] + I2[0])*2;  
		Q4[2] =  1                     - (I1[2] + I2[0])*2; 
	    
	    // These are used in the macro FILL_PPIPED
		int x0 = P1[0];
		int y0 = P1[1];
		int z0 = P1[2];

		O p1a[4], p1b[4], p2a[4], p2b[4], p3a[4], p3b[4], p4a[4], p4b[4];
		O p1a_r[4], p1b_r[4], p2a_r[4], p2b_r[4], p3a_r[4], p3b_r[4], p4a_r[4], p4b_r[4];

		// Do the 32 value lookups
		FILL_PPIPED(p1a[0], p1a[1], p1a[2], p1a[3],
					p1b[0], p1b[1], p1b[2], p1b[3], Q1);

		x0 = P1[0] + P2[0]; 
		y0 = P1[1] + P2[1]; 
		z0 = P1[2] + P2[2];

		FILL_PPIPED(p2a[0], p2a[1], p2a[2], p2a[3],
					p2b[0], p2b[1], p2b[2], p2b[3], Q2);

		x0 = P1[0] + P3[0]; 
		y0 = P1[1] + P3[1]; 
		z0 = P1[2] + P3[2];

		FILL_PPIPED(p3a[0], p3a[1], p3a[2], p3a[3],
					p3b[0], p3b[1], p3b[2], p3b[3], Q3);

		x0 = P1[0] + P4[0]; 
		y0 = P1[1] + P4[1]; 
		z0 = P1[2] + P4[2];
 
		FILL_PPIPED(p4a[0], p4a[1], p4a[2], p4a[3],
					p4b[0], p4b[1], p4b[2], p4b[3], Q4);

        // Permute the data points
        p1a_r[0] = I1[0] * p1a[0] + I1[1] * p1a[0] + I1[2] * p1a[0] +
                   I2[0] * p1a[0] + I2[1] * p1a[0] + I2[2] * p1a[0];
                   
        p1a_r[1] = I1[0] * p1a[1] + I1[1] * p1a[2] + I1[2] * p1a[2] +
                   I2[0] * p1a[1] + I2[1] * p1a[3] + I2[2] * p1a[3];
                   
        p1a_r[2] = I1[0] * p1a[2] + I1[1] * p1a[1] + I1[2] * p1a[3] + 
                   I2[0] * p1a[3] + I2[1] * p1a[1] + I2[2] * p1a[2];
                   
        p1a_r[3] = I1[0] * p1a[3] + I1[1] * p1a[3] + I1[2] * p1a[1] +
                   I2[0] * p1a[2] + I2[1] * p1a[2] + I2[2] * p1a[1];
                   
        // p1b
        p1b_r[0] = I1[0] * p1b[0] + I1[1] * p1b[1] + I1[2] * p1b[1] +
                   I2[0] * p1b[0] + I2[1] * p1b[2] + I2[2] * p1b[2];
                   
        p1b_r[1] = I1[0] * p1b[1] + I1[1] * p1b[0] + I1[2] * p1b[2] +
                   I2[0] * p1b[2] + I2[1] * p1b[0] + I2[2] * p1b[1];
                   
        p1b_r[2] = I1[0] * p1b[2] + I1[1] * p1b[2] + I1[2] * p1b[0] + 
                   I2[0] * p1b[1] + I2[1] * p1b[1] + I2[2] * p1b[0];
                   
        p1b_r[3] = I1[0] * p1b[3] + I1[1] * p1b[3] + I1[2] * p1b[3] +
                   I2[0] * p1b[3] + I2[1] * p1b[3] + I2[2] * p1b[3];
                   
        // p2a
        p2a_r[0] = I2[2] * p2a[0] + I2[1] * p2a[0] + I2[0] * p2a[0] +
                   I1[2] * p2a[0] + I1[1] * p2a[0] + I1[0] * p2a[0];
                   
        p2a_r[1] = I2[2] * p2a[1] + I2[1] * p2a[2] + I2[0] * p2a[2] +
                   I1[2] * p2a[1] + I1[1] * p2a[3] + I1[0] * p2a[3];
                   
        p2a_r[2] = I2[2] * p2a[2] + I2[1] * p2a[1] + I2[0] * p2a[3] + 
                   I1[2] * p2a[3] + I1[1] * p2a[1] + I1[0] * p2a[2];
                   
        p2a_r[3] = I2[2] * p2a[3] + I2[1] * p2a[3] + I2[0] * p2a[1] +
                   I1[2] * p2a[2] + I1[1] * p2a[2] + I1[0] * p2a[1];
                   
        // p2b
        p2b_r[0] = I2[2] * p2b[0] + I2[1] * p2b[1] + I2[0] * p2b[1] +
                   I1[2] * p2b[0] + I1[1] * p2b[2] + I1[0] * p2b[2];
                   
        p2b_r[1] = I2[2] * p2b[1] + I2[1] * p2b[0] + I2[0] * p2b[2] +
                   I1[2] * p2b[2] + I1[1] * p2b[0] + I1[0] * p2b[1];
                   
        p2b_r[2] = I2[2] * p2b[2] + I2[1] * p2b[2] + I2[0] * p2b[0] + 
                   I1[2] * p2b[1] + I1[1] * p2b[1] + I1[0] * p2b[0];
                   
        p2b_r[3] = I2[2] * p2b[3] + I2[1] * p2b[3] + I2[0] * p2b[3] +
                   I1[2] * p2b[3] + I1[1] * p2b[3] + I1[0] * p2b[3];
                   
        // p3a
        p3a_r[0] = I2[2] * p3a[0] + I2[1] * p3a[0] + I2[0] * p3a[0] +
                   I1[2] * p3a[0] + I1[1] * p3a[0] + I1[0] * p3a[0];
                   
        p3a_r[1] = I2[2] * p3a[1] + I2[1] * p3a[2] + I2[0] * p3a[2] +
                   I1[2] * p3a[1] + I1[1] * p3a[3] + I1[0] * p3a[3];
                   
        p3a_r[2] = I2[2] * p3a[2] + I2[1] * p3a[1] + I2[0] * p3a[3] + 
                   I1[2] * p3a[3] + I1[1] * p3a[1] + I1[0] * p3a[2];
                   
        p3a_r[3] = I2[2] * p3a[3] + I2[1] * p3a[3] + I2[0] * p3a[1] +
                   I1[2] * p3a[2] + I1[1] * p3a[2] + I1[0] * p3a[1];
                   
        // p3b
        p3b_r[0] = I2[2] * p3b[0] + I2[1] * p3b[1] + I2[0] * p3b[1] +
                   I1[2] * p3b[0] + I1[1] * p3b[2] + I1[0] * p3b[2];
                   
        p3b_r[1] = I2[2] * p3b[1] + I2[1] * p3b[0] + I2[0] * p3b[2] +
                   I1[2] * p3b[2] + I1[1] * p3b[0] + I1[0] * p3b[1];
                   
        p3b_r[2] = I2[2] * p3b[2] + I2[1] * p3b[2] + I2[0] * p3b[0] + 
                   I1[2] * p3b[1] + I1[1] * p3b[1] + I1[0] * p3b[0];
                   
        p3b_r[3] = I2[2] * p3b[3] + I2[1] * p3b[3] + I2[0] * p3b[3] +
                   I1[2] * p3b[3] + I1[1] * p3b[3] + I1[0] * p3b[3];
        
         // p4a
        p4a_r[0] = I1[0] * p4a[0] + I1[1] * p4a[0] + I1[2] * p4a[0] +
                   I2[0] * p4a[0] + I2[1] * p4a[0] + I2[2] * p4a[0];
                   
        p4a_r[1] = I1[0] * p4a[1] + I1[1] * p4a[2] + I1[2] * p4a[2] +
                   I2[0] * p4a[1] + I2[1] * p4a[3] + I2[2] * p4a[3];
                   
        p4a_r[2] = I1[0] * p4a[2] + I1[1] * p4a[1] + I1[2] * p4a[3] + 
                   I2[0] * p4a[3] + I2[1] * p4a[1] + I2[2] * p4a[2];
                   
        p4a_r[3] = I1[0] * p4a[3] + I1[1] * p4a[3] + I1[2] * p4a[1] +
                   I2[0] * p4a[2] + I2[1] * p4a[2] + I2[2] * p4a[1];
                   
        // p4b
        p4b_r[0] = I1[0] * p4b[0] + I1[1] * p4b[1] + I1[2] * p4b[1] +
                   I2[0] * p4b[0] + I2[1] * p4b[2] + I2[2] * p4b[2];
                   
        p4b_r[1] = I1[0] * p4b[1] + I1[1] * p4b[0] + I1[2] * p4b[2] +
                   I2[0] * p4b[2] + I2[1] * p4b[0] + I2[2] * p4b[1];
                   
        p4b_r[2] = I1[0] * p4b[2] + I1[1] * p4b[2] + I1[2] * p4b[0] + 
                   I2[0] * p4b[1] + I2[1] * p4b[1] + I2[2] * p4b[0];
                   
        p4b_r[3] = I1[0] * p4b[3] + I1[1] * p4b[3] + I1[2] * p4b[3] +
                   I2[0] * p4b[3] + I2[1] * p4b[3] + I2[2] * p4b[3];

		float a = maxParameter - 1.0, b = midParameter - 1.0, g = minParameter - 1.0;
		CONV_PPIPED(r, p1a_r, p1b_r, a, b, g);

		a = -minParameter; b = maxParameter-minParameter-1.0; g = midParameter-minParameter-1.0;
		CONV_PPIPED(r, p2a_r, p2b_r, a, b, g);

		a = -maxParameter+midParameter; b = -maxParameter+minParameter; g = -maxParameter;
		CONV_PPIPED(r, p3a_r, p3b_r, a, b, g);
        
		a = -midParameter+minParameter; b = -midParameter; g = maxParameter-midParameter-1.0;
		CONV_PPIPED(r, p4a_r, p4b_r, a, b, g);

		return r;
#undef RHO_FAST
#undef RHO22
#undef FILL_PPIPED
#undef FILL_PPIPED_NORMAL
#undef CONV_PPIPED
	}

	static const vector3<O> convolutionSumNormal(const vector3<I> &p, shift_invariant_space3<quintic_box, O, I> *L) {
//
#define RHO_FAST(a, b, g, a2, a3) \
    -a3*(g*b-0.5*a*(g+b)+0.3*a2);

#define RHO22(A, B, G) ( -((A)*(A))*(A) \
    * ((1./6.)*(G)*(B) - (1./12.) *(A) \
    * ((G)+(B))+0.05 * ((A)*(A))))

#define FILL_PPIPED(p0, p1, p2, p3, p4, p5, p6, p7, Q) { \
	p0 = L->GN(x0,        y0,        z0      ); \
	p1 = L->GN(x0 - Q[0], y0 + Q[1], z0 + Q[2]); \
	p2 = L->GN(x0 + Q[0], y0 - Q[1], z0 + Q[2]); \
	p3 = L->GN(x0 + Q[0], y0 + Q[1], z0 - Q[2]); \
	p4 = L->GN(x0 + 2*Q[0], y0,          z0         ); \
	p5 = L->GN(x0,          y0 + 2*Q[1], z0         ); \
	p6 = L->GN(x0,          y0,          z0 + 2*Q[2]); \
	p7 = L->GN(x0 + Q[0],   y0 + Q[1],   z0 + Q[2]   ); \
	}

#define CONV_PPIPED(value, pa, pb, alpha, beta, gamma) { \
    vector3<O> p123 = pa[1] + pa[2] + pa[3]; \
    vector3<O> p0 = pa[0] * 4.0; \
    O alpha2 = alpha * alpha; \
    O alpha3 = alpha2 * alpha * (1./6.); \
    O base_rho = RHO_FAST(alpha, beta, gamma, alpha2, alpha3); \
    O base_rho2 = base_rho - 0.5 * alpha3 * alpha; \
    value = value + (p0*(-2.5)+(p123)*4.0-(pb[0]+pb[1]+pb[2])*2.0 + pb[3]) * base_rho; \
    alpha2 = base_rho2 + alpha3 * beta; \
    base_rho2 += alpha3 * gamma; \
    value = value + (p0-(p123-pa[1])*2.0+pb[0]) * (alpha2) \
    + (p0-(p123-pa[2])*2.0+pb[1]) * (base_rho2) \
    + (p0*(-.5)+pa[3]) \
    * (base_rho2 + alpha2 - alpha3 - base_rho); \
    alpha -= 1.0; alpha2 = beta * beta; \
    alpha3 = (1./6.) * alpha2 * beta; \
    base_rho = RHO_FAST(beta, gamma, alpha, alpha2, alpha3); \
    value = value + (p0-(p123-pa[3])*2.0+pb[2]) * (base_rho); \
    p0 = p0*(-.5); \
    value = value + (p0+pa[2]) * (base_rho + alpha3 * (alpha - .5*beta)) \
    + (p0 + pa[1]) * RHO22(gamma, alpha, beta-1.0) + (p0 * (-0.5) ) \
    * RHO22(alpha, beta-1.0, gamma-1.0); \
	}

		int P1[3],P2[3],P3[3],P4[3];
		int Q1[3],Q2[3],Q3[3],Q4[3];
		int I1[3],I2[3];

		vector3<O> r = vector3<O>(0,0,0);
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

		// We have six cases, the above sort decides which we 
		// use, later on we permute according to this
		I1[0] = (i == 7); I2[0] = (i == 3);
		I1[1] = (i == 5); I2[1] = (i == 2);
		I1[2] = (i == 4); I2[2] = (i == 0);

		// The first two points come for free
		// diagonal
		P2[0] = 1; Q1[0] = -1;
		P2[1] = 1; Q1[1] = -1;
		P2[2] = 1; Q1[2] = -1;
	
		// The rest are determinied by the sort
		Q3[0] = P3[0] =  1 - (I1[2] + I2[2])*2; 
		Q3[1] = P3[1] =  1                     - (I2[0] + I2[1])*2; 
		Q3[2] = P3[2] = -1 + (I1[2] + I2[2])*2 + (I2[0] + I2[1])*2;  

		P4[0] =  2 - (I1[1] + I1[2])*2 - (I2[1] + I2[2])*2; 
		P4[1] =  0 + (I1[1] + I1[2])*2;  
		P4[2] =  0                     + (I2[1] + I2[2])*2; 

		Q2[0] = -1 + (I1[1] + I1[2])*2 + (I2[1] + I2[2])*2; 
		Q2[1] =  1 - (I1[1] + I1[2])*2;  
		Q2[2] =  1                     - (I2[1] + I2[2])*2; 

		Q4[0] =  1 - (I1[1] + I2[1])*2 ; 
		Q4[1] = -1 + (I1[1] + I2[1])*2 + (I1[2] + I2[0])*2;  
		Q4[2] =  1                     - (I1[2] + I2[0])*2; 
	    
	    // These are used in the macro FILL_PPIPED
		int x0 = P1[0];
		int y0 = P1[1];
		int z0 = P1[2];

		vector3<O> p1a[4], p1b[4], p2a[4], p2b[4], p3a[4], p3b[4], p4a[4], p4b[4];
		vector3<O> p1a_r[4], p1b_r[4], p2a_r[4], p2b_r[4], p3a_r[4], p3b_r[4], p4a_r[4], p4b_r[4];

		// Do the 32 value lookups
		FILL_PPIPED(p1a[0], p1a[1], p1a[2], p1a[3],
					p1b[0], p1b[1], p1b[2], p1b[3], Q1);

		x0 = P1[0] + P2[0]; 
		y0 = P1[1] + P2[1]; 
		z0 = P1[2] + P2[2];

		FILL_PPIPED(p2a[0], p2a[1], p2a[2], p2a[3],
					p2b[0], p2b[1], p2b[2], p2b[3], Q2);

		x0 = P1[0] + P3[0]; 
		y0 = P1[1] + P3[1]; 
		z0 = P1[2] + P3[2];

		FILL_PPIPED(p3a[0], p3a[1], p3a[2], p3a[3],
					p3b[0], p3b[1], p3b[2], p3b[3], Q3);

		x0 = P1[0] + P4[0]; 
		y0 = P1[1] + P4[1]; 
		z0 = P1[2] + P4[2];
 
		FILL_PPIPED(p4a[0], p4a[1], p4a[2], p4a[3],
					p4b[0], p4b[1], p4b[2], p4b[3], Q4);
        // Permute the data points
        p1a_r[0] = p1a[0] * I1[0] + p1a[0] * I1[1] + p1a[0] * I1[2] +
                   p1a[0] * I2[0] + p1a[0] * I2[1] + p1a[0] * I2[2];
                   
        p1a_r[1] = p1a[1] * I1[0] + p1a[2] * I1[1] + p1a[2] * I1[2] +
                   p1a[1] * I2[0] + p1a[3] * I2[1] + p1a[3] * I2[2];
                   
        p1a_r[2] = p1a[2] * I1[0] + p1a[1] * I1[1] + p1a[3] * I1[2] + 
                   p1a[3] * I2[0] + p1a[1] * I2[1] + p1a[2] * I2[2];
                   
        p1a_r[3] = p1a[3] * I1[0] + p1a[3] * I1[1] + p1a[1] * I1[2] +
                   p1a[2] * I2[0] + p1a[2] * I2[1] + p1a[1] * I2[2];
                   
        // p1b
        p1b_r[0] = p1b[0] * I1[0] + p1b[1] * I1[1] + p1b[1] * I1[2] +
                   p1b[0] * I2[0] + p1b[2] * I2[1] + p1b[2] * I2[2];
                   
        p1b_r[1] = p1b[1] * I1[0] + p1b[0] * I1[1] + p1b[2] * I1[2] +
                   p1b[2] * I2[0] + p1b[0] * I2[1] + p1b[1] * I2[2];
                   
        p1b_r[2] = p1b[2] * I1[0] + p1b[2] * I1[1] + p1b[0] * I1[2] + 
                   p1b[1] * I2[0] + p1b[1] * I2[1] + p1b[0] * I2[2];
                   
        p1b_r[3] = p1b[3] * I1[0] + p1b[3] * I1[1] + p1b[3] * I1[2] +
                   p1b[3] * I2[0] + p1b[3] * I2[1] + p1b[3] * I2[2];
                   
        // p2a
        p2a_r[0] = p2a[0] * I2[2] + p2a[0] * I2[1] + p2a[0] * I2[0] +
                   p2a[0] * I1[2] + p2a[0] * I1[1] + p2a[0] * I1[0];
                   
        p2a_r[1] = p2a[1] * I2[2] + p2a[2] * I2[1] + p2a[2] * I2[0] +
                   p2a[1] * I1[2] + p2a[3] * I1[1] + p2a[3] * I1[0];
                   
        p2a_r[2] = p2a[2] * I2[2] + p2a[1] * I2[1] + p2a[3] * I2[0] + 
                   p2a[3] * I1[2] + p2a[1] * I1[1] + p2a[2] * I1[0];
                   
        p2a_r[3] = p2a[3] * I2[2] + p2a[3] * I2[1] + p2a[1] * I2[0] +
                   p2a[2] * I1[2] + p2a[2] * I1[1] + p2a[1] * I1[0];
                   
        // p2b
        p2b_r[0] = p2b[0] * I2[2] + p2b[1] * I2[1] + p2b[1] * I2[0] +
                   p2b[0] * I1[2] + p2b[2] * I1[1] + p2b[2] * I1[0];
                   
        p2b_r[1] = p2b[1] * I2[2] + p2b[0] * I2[1] + p2b[2] * I2[0] +
                   p2b[2] * I1[2] + p2b[0] * I1[1] + p2b[1] * I1[0];
                   
        p2b_r[2] = p2b[2] * I2[2] + p2b[2] * I2[1] + p2b[0] * I2[0] + 
                   p2b[1] * I1[2] + p2b[1] * I1[1] + p2b[0] * I1[0];
                   
        p2b_r[3] = p2b[3] * I2[2] + p2b[3] * I2[1] + p2b[3] * I2[0] +
                   p2b[3] * I1[2] + p2b[3] * I1[1] + p2b[3] * I1[0];
                   
        // p3a
        p3a_r[0] = p3a[0] * I2[2] + p3a[0] * I2[1] + p3a[0] * I2[0] +
                   p3a[0] * I1[2] + p3a[0] * I1[1] + p3a[0] * I1[0];
                   
        p3a_r[1] = p3a[1] * I2[2] + p3a[2] * I2[1] + p3a[2] * I2[0] +
                   p3a[1] * I1[2] + p3a[3] * I1[1] + p3a[3] * I1[0];
                   
        p3a_r[2] = p3a[2] * I2[2] + p3a[1] * I2[1] + p3a[3] * I2[0] + 
                   p3a[3] * I1[2] + p3a[1] * I1[1] + p3a[2] * I1[0];
                   
        p3a_r[3] = p3a[3] * I2[2] + p3a[3] * I2[1] + p3a[1] * I2[0] +
                   p3a[2] * I1[2] + p3a[2] * I1[1] + p3a[1] * I1[0];
                   
        // p3b
        p3b_r[0] = p3b[0] * I2[2] + p3b[1] * I2[1] + p3b[1] * I2[0] +
                   p3b[0] * I1[2] + p3b[2] * I1[1] + p3b[2] * I1[0];
                   
        p3b_r[1] = p3b[1] * I2[2] + p3b[0] * I2[1] + p3b[2] * I2[0] +
                   p3b[2] * I1[2] + p3b[0] * I1[1] + p3b[1] * I1[0];
                   
        p3b_r[2] = p3b[2] * I2[2] + p3b[2] * I2[1] + p3b[0] * I2[0] + 
                   p3b[1] * I1[2] + p3b[1] * I1[1] + p3b[0] * I1[0];
                   
        p3b_r[3] = p3b[3] * I2[2] + p3b[3] * I2[1] + p3b[3] * I2[0] +
                   p3b[3] * I1[2] + p3b[3] * I1[1] + p3b[3] * I1[0];
        
         // p4a
        p4a_r[0] = p4a[0] * I1[0] + p4a[0] * I1[1] + p4a[0] * I1[2] +
                   p4a[0] * I2[0] + p4a[0] * I2[1] + p4a[0] * I2[2];
                   
        p4a_r[1] = p4a[1] * I1[0] + p4a[2] * I1[1] + p4a[2] * I1[2] +
                   p4a[1] * I2[0] + p4a[3] * I2[1] + p4a[3] * I2[2];
                   
        p4a_r[2] = p4a[2] * I1[0] + p4a[1] * I1[1] + p4a[3] * I1[2] + 
                   p4a[3] * I2[0] + p4a[1] * I2[1] + p4a[2] * I2[2];
                   
        p4a_r[3] = p4a[3] * I1[0] + p4a[3] * I1[1] + p4a[1] * I1[2] +
                   p4a[2] * I2[0] + p4a[2] * I2[1] + p4a[1] * I2[2];
                   
        // p4b
        p4b_r[0] = p4b[0] * I1[0] + p4b[1] * I1[1] + p4b[1] * I1[2] +
                   p4b[0] * I2[0] + p4b[2] * I2[1] + p4b[2] * I2[2];
                   
        p4b_r[1] = p4b[1] * I1[0] + p4b[0] * I1[1] + p4b[2] * I1[2] +
                   p4b[2] * I2[0] + p4b[0] * I2[1] + p4b[1] * I2[2];
                   
        p4b_r[2] = p4b[2] * I1[0] + p4b[2] * I1[1] + p4b[0] * I1[2] + 
                   p4b[1] * I2[0] + p4b[1] * I2[1] + p4b[0] * I2[2];
                   
        p4b_r[3] = p4b[3] * I1[0] + p4b[3] * I1[1] + p4b[3] * I1[2] +
                   p4b[3] * I2[0] + p4b[3] * I2[1] + p4b[3] * I2[2];

		float a = maxParameter - 1.0, b = midParameter - 1.0, g = minParameter - 1.0;
		CONV_PPIPED(r, p1a_r, p1b_r, a, b, g);

		a = -minParameter; b = maxParameter-minParameter-1.0; g = midParameter-minParameter-1.0;
		CONV_PPIPED(r, p2a_r, p2b_r, a, b, g);

		a = -maxParameter+midParameter; b = -maxParameter+minParameter; g = -maxParameter;
		CONV_PPIPED(r, p3a_r, p3b_r, a, b, g);
        
		a = -midParameter+minParameter; b = -midParameter; g = maxParameter-midParameter-1.0;
		CONV_PPIPED(r, p4a_r, p4b_r, a, b, g);

		return r;
#undef RHO_FAST
#undef RHO22
#undef FILL_PPIPED
#undef FILL_PPIPED_NORMAL
#undef CONV_PPIPED
	}

	//Filters	
	static O interpolationFilter(const I &u, const I &v, const I &w, const I &h) {
		return 1./((1./15.) * (6. + 6.*cos(u)*cos(v)*cos(w) + cos(2.*u) + cos(2.*v) + cos(2.*w)));
	}

	static O poissonFilter(const I &u, const I &v, const I &w, const I &h) {
		O invh = 1./(h*h);
		return 1./(invh*(1./6.) * (-15. + 16.*cos(u)*cos(v)*cos(w) - cos(2.*u)*cos(2.*v)*cos(2.*w)));
	}

#define bcc_odd_divFilter(O,L,x,y,z) \
	 [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){ \
			double iscl = -1./(((L*)x)->getScale()*12.); \
			return 0.; \
		} \

	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		divergenceFilter(L *x, L *y, L *z){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			printf("0x%x, 0x%x, 0x%x\n", x, y, z);
//			double iscl = -1./(x->getScale()*12.);

/*			double dx = 
				 -1.*x->GV(i-2, j+2, k+2) + 
				  8.*x->GV(i-1, j+1, k+1) +
				 -8.*x->GV(i+1, j-1, k-1) + 
					 x->GV(i+2, j-2, k-2);
			double dy = 
				 -1.*y->GV(i+2, j-2, k+2) + 
				  8.*y->GV(i+1, j-1, k+1) +
				 -8.*y->GV(i-1, j+1, k-1) + 
					 y->GV(i-2, j+2, k-2);
					
			double dz = 
				 -1.*z->GV(i+2, j+2, k-2) + 
				  8.*z->GV(i+1, j+1, k-1) +
				 -8.*z->GV(i-1, j-1, k+1) + 
					 z->GV(i-2, j-2, k+2);
*/
			return 0;// iscl*(dx+dy+dz);
		};
	}


	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		shiftedDivergenceFilter(L *x, L *y, L *z){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			double iscl = -1./(x->getScale()*24.);

			double dx = 
				 -1.*x->GV(i-2, j+2, k+2) + 
				  27.*x->GV(i-1, j+1, k+1) +
				 -27.*x->GV(i, j, k) + 
					 x->GV(i+1, j-1, k-1);

			double dy = 
				 -1.*y->GV(i+2, j-2, k+2) + 
				  27.*y->GV(i+1, j-1, k+1) +
				 -27.*y->GV(i, j, k) + 
					 y->GV(i-1, j+1, k-1);
					
			double dz = 
				 -1.*z->GV(i+2, j+2, k-2) + 
				  27.*z->GV(i+1, j+1, k-1) +
				 -27.*z->GV(i, j, k) + 
					 z->GV(i-1, j-1, k+1);

			return iscl*(dx+dy+dz);
		};
	}
};

};
#endif