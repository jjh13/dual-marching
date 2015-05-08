#include <iostream>
#include <tclap/CmdLine.h>

// SISL Main Include
#include <sisl/sisl.hpp>

// odd cartesian function spaces
#include <sisl/lattice/cartesian_odd.hpp>
#include <sisl/basis/tp3cubic.hpp>

// odd BCC function spaces
#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/basis/quintic.hpp>
#include <sisl/basis/linear_rdod.hpp>

// utility functions
#include <poisson/pointset.hpp>
#include <sisl/utility/isosurface.hpp>
#include <sisl/utility/dualbcc.hpp>
#include <sisl/utility/dualfcc.hpp>
#include <sisl/utility/dualcc.hpp>
#include <sisl/utility/scattered.hpp>
#include <sisl/utility/ply_writer.hpp>

#include <Eigen/Dense>

#include <tuple>
#include <cmath>

#define VESION_STRING "0.1"

using namespace sisl;
using namespace std;
using namespace TCLAP;

template <class T>
class MarschnerLobb {
public:
	double a, f;

	MarschnerLobb(double Fm, double alpha) : f(Fm), a(alpha){ }
	double rho(double r){
		return cos(2*M_PI*f*cos(r*M_PI/2.));
	}	
	double evaluate(const double &xx, const double &yy, const double &zz) {
		double x = (xx-0.5)*2., y = (yy-0.5)*2., z = (zz-0.5)*2.;

		double r = rho(sqrt(x*x + y*y));


		double ret = 1. - sin(M_PI*z/2.) + a * (1. + r);

		return ret/(2. + 2.*a);

	}

	double evaluate(const vector3<T> &p) {return this->evaluate(p.i, p.j, p.k);}
	double div(const vector3<T> &p) { return this->div(p.i, p.j, p.k); }
	double div(const double &x, const double &y, const double &z) {return 0;}
	double laplacian(const vector3<T> &p){return this->laplacian(p.i, p.j, p.k);}
	double laplacian(const double &x, const double &y, const double &z) {return 0;}
};

template <class T>
class HamFunction {
public:
	HamFunction(){}
	double evaluate(const double &x, const double &y, const double &z) {
		// double xx = (x-0.5), yy = (y-0.5), zz = (z-0.5);
		// return xx*xx + yy*yy + zz*zz - 0.25*0.25;

		if( (x >= 0.25 && x <= 0.75) && (y >= 0.25 && y <= 0.75) && (z >= 0.25 && z <= 0.75)) {
			return -1;
		}
		return 1.;
		// return -(sin(0.3141592654e1 * x) * sin(0.3141592654e1 * y) * sin(0.3141592654e1 * z) * (sqrt(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1)) - 0.2e1 * cos(0.8e1 * 0.3141592654e1 * (0.9e1 * z - 0.45e1) * pow(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1), -0.1e1 / 0.2e1)) - 0.2e1));
	}

	double evaluate(const vector3<T> &p) {
		return this->evaluate(p.i, p.j, p.k);
	}
	
	double div(const vector3<T> &p) { return this->div(p.i, p.j, p.k); }
	double div(const double &x, const double &y, const double &z) {
		return 0;
	}
	
	double laplacian(const vector3<T> &p){return this->laplacian(p.i, p.j, p.k);}
	double laplacian(const double &x, const double &y, const double &z) {
		return 0;
	}
};

int main(int argc, char *argv[])
{
	try {
		CmdLine cmd("Dual BCC marching", ' ', VESION_STRING);
		sisl::utility::marchingCubes<double> mc;
		sisl::utility::dualbcc_isosurface<double> dbcc;
		sisl::utility::dualfcc_isosurface<double> dfcc;
		sisl::utility::dualcc_isosurface<double> dcc;

		ValueArg<std::string> outputArg("o", "output", "Output mesh/implicit file",true,"out.ply","filename");

		//cmd.add(outputArg);
		cmd.parse(argc, argv);
//		std::string output = outputArg.getValue();

		typedef linear_bcc_box<double, double> GFType;
		typedef bcc_odd<linear_bcc_box<double, double>, double, double> LATType;

		HamFunction<double> f;
		// MarschnerLobb<double> f(6, 0.25);

		LATType *lat = new LATType(1./double(2.*50));

		lat->forEachLatticeSite([&](const int &i, const int &j, const int &k) {
			vector3<double> p = lat->getSitePosition(i,j,k);
			return f.evaluate(p);
		});

		dbcc.contour<LATType, double, double>(
			lat, 0.5, 1./(2.*25),
			vector3<double>(0.1,0.1,0.1), 
			vector3<double>(0.9,0.9,0.9)
		);
		dfcc.contour<LATType, double, double>(
			lat, 0.5, 1./(2.*25),
			vector3<double>(0.1,0.1,0.1), 
			vector3<double>(0.9,0.9,0.9)
		);
		dcc.contour<LATType, double, double>(
			lat, 0.5, 1./(2.*25),
			vector3<double>(0.1,0.1,0.1), 
			vector3<double>(0.9,0.9,0.9)
		);

		dbcc.writeSurface("dualBCCTest.ply");
		dfcc.writeSurface("dualFCCTest.ply");
		dcc.writeSurface("dualCCTest.ply");
		// mc.writeSurface("marchingCubesTest.ply");

	}catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
	}catch (char const* e) {
		cerr << e << endl; 
	}
}


