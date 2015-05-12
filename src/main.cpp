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
	double a, fm;

	MarschnerLobb(double Fm, double alpha) : fm(Fm), a(alpha){ }
	double rho(double r){
		return cos(2*M_PI*fm*cos(r*M_PI/2.));
	}	
	double f(const double &xx, const double &yy, const double &zz) {
		double x = (2.*xx-1.), y = (2.*yy-1.), z = (2.*zz-1.);
		double r = rho(sqrt(x*x + y*y));
		double ret = 1. - sin(M_PI*z/2.) + a * (1. + r);
		return ret/(2. + 2.*a);

	}

	double f(const vector3<T> &p) {return this->f(p.i, p.j, p.k);}
	vector3<double> grad_f(const vector3<T> &p) { return this->grad_f(p.i, p.j, p.k); }
	vector3<double> grad_f(const double &x, const double &y, const double &z) {
		double xx = 2.*x - 1.;
		double yy = 2.*y - 1.;
		double zz = 2.*z - 1.;
		return vector3<double> (
			M_PI*M_PI*a*fm*(xx)*sin(2.*M_PI*fm*cos(0.5*M_PI*sqrt(xx*xx + yy*yy))) *sin(0.5*M_PI*sqrt(xx*xx + yy*yy))/(sqrt(xx*xx + yy*yy)*(a + 1.)),
			M_PI*M_PI*a*fm*(yy)*sin(2.*M_PI*fm*cos(0.5*M_PI*sqrt(xx*xx + yy*yy))) *sin(0.5*M_PI*sqrt(xx*xx + yy*yy))/(sqrt(xx*xx + yy*yy)*(a + 1.)),
			-0.5*M_PI*cos(0.5*M_PI*zz)/(a + 1.)
		);
	}
};

template <class T>
class HamFunction {
public:
	HamFunction(){}
	double f(const double &x, const double &y, const double &z) {
		return -(sin(0.3141592654e1 * x) * sin(0.3141592654e1 * y) * sin(0.3141592654e1 * z) * (sqrt(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1)) - 0.2e1 * cos(0.8e1 * 0.3141592654e1 * (0.9e1 * z - 0.45e1) * pow(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1), -0.1e1 / 0.2e1)) - 0.2e1));
	}

	double f(const vector3<T> &p) {
		return this->f(p.i, p.j, p.k);
	}
	
	vector3<double> grad_f(const vector3<T> &p) { return this->grad_f(p.i, p.j, p.k); }
	vector3<double> grad_f(const double &x, const double &y, const double &z) {
		return vector3<double> (
			2.*(x-0.5),
			2.*(y-0.5),
			2.*(z-0.5)
		);
	}
};


template<class T>
class SphereFunction{
public:
	SphereFunction(){}
	double f(const double &x, const double &y, const double &z) {
		double xx = (x-0.5), yy = (y-0.5), zz = (z-0.5);
		return xx*xx + yy*yy + zz*zz - 0.25*0.25;
	}

	double f(const vector3<T> &p) {return this->f(p.i, p.j, p.k);}

	vector3<double> grad_f(const vector3<T> &p) { return this->grad_f(p.i, p.j, p.k); }
	vector3<double> grad_f(const double &x, const double &y, const double &z) {
		return vector3<double> (
			2*(x-0.5),
			2*(y-0.5),
			2*(z-0.5)
		);
	}
};

int main(int argc, char *argv[])
{
	try {
		CmdLine cmd("Dual marching cubes for CC/BCC/FCC lattices", ' ', VESION_STRING);
		sisl::utility::marchingCubes<double> mc;
		sisl::utility::dualbcc_isosurface<double> dbcc;
		sisl::utility::dualfcc_isosurface<double> dfcc;
		sisl::utility::dualcc_isosurface<double> dcc;

		ValueArg<std::string> outputArg("o", "output", "Output mesh name", true,"output", "filename");
		ValueArg<std::string> testFunction("t", "test_function", "test function", false, "lobner", "test function");
		ValueArg<double> isoValue("i", "iso_value", "Isovalue for contour", true, 0, "iso-value");
		ValueArg<double> gridGranularity("s", "grid_granularity", "Grid grid granularity",true, 0.25, "dh");

		cmd.add(testFunction);
		cmd.add(outputArg);
		cmd.add(isoValue);
		cmd.add(gridGranularity);

		cmd.parse(argc, argv);
		std::string output = outputArg.getValue();
		std::string function = testFunction.getValue();

		double levelset = isoValue.getValue();

		HamFunction<double> hf;
		SphereFunction<double> sf;
		MarschnerLobb<double> mf(6, 0.25);

		if(function == std::string("lobb")){
			dbcc.contour<MarschnerLobb<double>, double, double>(
				&mf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
			dfcc.contour<MarschnerLobb<double>, double, double>(
				&mf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
			dcc.contour<MarschnerLobb<double>, double, double>(
				&mf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
		}else if(function == std::string("ham")){
			dbcc.contour<HamFunction<double>, double, double>(
				&hf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
			dfcc.contour<HamFunction<double>, double, double>(
				&hf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
			dcc.contour<HamFunction<double>, double, double>(
				&hf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
		}else if(function == std::string("sphere")){
			dbcc.contour<SphereFunction<double> , double, double>(
				&sf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
			dfcc.contour<SphereFunction<double> , double, double>(
				&sf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
			dcc.contour<SphereFunction<double> , double, double>(
				&sf, levelset, 1./(2.*100),
				vector3<double>(0,0,0), 
				vector3<double>(1,1,1)
			);
		}


		dbcc.writeSurface(output + std::string(".bcc.ply"));
		dfcc.writeSurface(output + std::string(".fcc.ply"));
		dcc.writeSurface(output + std::string(".cc.ply"));

	}catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
	}catch (char const* e) {
		cerr << e << endl; 
	}
}


