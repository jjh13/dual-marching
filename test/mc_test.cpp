#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>

#include <iostream>

#include <sisl/utility/isosurface.hpp>

using namespace sisl;

class testSurface {
public:
	testSurface() {}
	virtual double f(const double &x, const double &y, const double &z) const {
    	return (sin(0.3141592654e1 * x) * sin(0.3141592654e1 * y) * sin(0.3141592654e1 * z) * (sqrt(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1)) - 0.2e1 * cos(0.8e1 * 0.3141592654e1 * (0.9e1 * z - 0.45e1) * pow(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1), -0.1e1 / 0.2e1)) - 0.2e1));
	}

	virtual double f(const vector3<double> &p) const {
		return f(p.i, p.j, p.k);
	}

	virtual vector3<double> grad_f(const vector3<double> &p) const {
		return vector3<double>(0,0,0);
	}

	virtual bool saveAsImplicit(const std::string &filename, bool writeGrad = false, double levelset = 0.) const {return false;}

	virtual void prefilter(bool interp = true, bool poisson = true, bool autocor = true) {}
	virtual sisl::vector3<double> getPosition(const int &x, const int &y, const int &z) const {
			double dh  = 1/2.;
			return vector3<double>(double(x)*dh,double(y)*dh,double(z)*dh); 
		};
	virtual double GV(const int &x, const int &y, const int &z) const{return 0;};
	virtual void SV(const int &x, const int &y, const int &z, const double &value) {};
	virtual int lIndex(const int &x, const int &y, const int &z) const {return -1;}; 
	virtual void *getArray(){return NULL;}
	virtual void setShift(const sisl::vector3<double> &){};
	virtual void getLatticeSite(const sisl::vector3<double> &, int &, int &, int &) const{};
};

TEST_CASE("Marching Cubes Test","mctest"){
	sisl::utility::marchingCubes<double> mc;
	testSurface t;
	try{
		mc.marchLattice<testSurface, double, double>(
			&t, 
			NULL, 
			NULL, 
			NULL, 
			0.2, 
			0.005, 
			vector3<double>(0,0,0), 
			vector3<double>(1,1,1));

		REQUIRE(mc.writeSurface("marchingCubesTest.ply"));
	}catch(char const*  c){
		std::cout << c << std::endl;
	}
}