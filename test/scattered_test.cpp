#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>

#include <sisl/lattice/cartesian_odd.hpp>
#include <sisl/basis/tp3cubic.hpp>

#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/basis/quintic.hpp>

#include <sisl/utility/scattered.hpp>
#include <sisl/utility/isosurface.hpp>

#include <iostream>

#define N_MC_SAMPLES  100000
#define N_VAR_SAMPLES 300000

using namespace sisl;

template <class T>
class HamFunction {
public:
	HamFunction(){}	
	double evaluate(const double &x, const double &y, const double &z) {
		return (sin(0.3141592654e1 * x) * sin(0.3141592654e1 * y) * sin(0.3141592654e1 * z) * (sqrt(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1)) - 0.2e1 * cos(0.8e1 * 0.3141592654e1 * (0.9e1 * z - 0.45e1) * pow(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1), -0.1e1 / 0.2e1)) - 0.2e1));
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

TEST_CASE("BCC Scattered Splatting Reconstruction", "bccsplatr"){
	typedef sisl::quintic_box<double,double> gType;
	typedef bcc_odd<quintic_box<double,double>, double, double> lType;

	utility::marchingCubes<double> mc;
	HamFunction<double> f;
	bcc_odd<quintic_box<double,double>, double, double> lat(1./double(2.*100));

//	viennacl::ocl::set_context_device_type(0, viennacl::ocl::cpu_tag());
//	printf("OpenCL Info:\n%s \n", viennacl::ocl::current_device().info().c_str());

	// Make a pointset
	std::vector<vector3<double>> pts(N_VAR_SAMPLES*10);
	std::vector<double> smpls(N_VAR_SAMPLES*10);

	#pragma omp parallel for
	for(int i = 0; i < N_VAR_SAMPLES*10; i++){
		vector3<double>  p(	((double)rand()/(double)RAND_MAX),
							((double)rand()/(double)RAND_MAX),
							((double)rand()/(double)RAND_MAX));
		pts[i] = p;
		smpls[i] = f.evaluate(p);
	}

	utility::splatSamplesOntoLattice<lType, gType, double, double>(
				pts,
				smpls,
				&lat,
				vector3<double>(0,0,0));

	lat.frequencyFilter(quintic_box<double,double>::interpolationFilter);

	mc.marchLattice<lType, double, double>(
		&lat, 
		NULL, 
		NULL, 
		NULL, 
		0.5, 
		0.01, 
		vector3<double>(0,0,0), 
		vector3<double>(1,1,1));
	
	REQUIRE(mc.writeSurface("scattered_splat_bcc_test.ply"));
}

TEST_CASE("Cartesian Scattered Data Reconstruction", "bccscatr"){
	typedef sisl::tp3cubic<double,double> gType;
	typedef cartesian_odd<tp3cubic<double,double>, double, double> lType;

	utility::marchingCubes<double> mc;
	HamFunction<double> f;
	lType lat(1./double(127));

//	viennacl::ocl::set_context_device_type(0, viennacl::ocl::cpu_tag());
//	printf("OpenCL Info:\n%s \n", viennacl::ocl::current_device().info().c_str());

	// Make a pointset
	std::vector<vector3<double>> pts(N_VAR_SAMPLES);
	std::vector<double> smpls(N_VAR_SAMPLES);

	#pragma omp parallel for
	for(int i = 0; i < N_VAR_SAMPLES; i++){
		vector3<double>  p(	((double)rand()/(double)RAND_MAX),
							((double)rand()/(double)RAND_MAX),
							((double)rand()/(double)RAND_MAX));
		pts[i] = p;
		smpls[i] = f.evaluate(p);
	}
	std::cout << sisl::utility::_cc_optimization_kernel_dbl << std :: endl;

	utility::fitPoints<lType, gType, double, double> (
				pts,
				smpls,
				gType::getBeppoLevi2Norm(),
				&lat,
				vector3<double>(0,0,0)
			);
	printf("About to march\n");

	mc.marchLattice<lType, double, double>(
		&lat, 
		NULL, 
		NULL, 
		NULL, 
		0.5, 
		0.01, 
		vector3<double>(0,0,0), 
		vector3<double>(1,1,1));
	REQUIRE(mc.writeSurface("scattered_cc_test.ply"));
}

TEST_CASE("BCC Scattered Data Reconstruction", "bccscatr"){
	typedef sisl::quintic_box<double,double> gType;
	typedef bcc_odd<quintic_box<double,double>, double, double> lType;

	utility::marchingCubes<double> mc;
	HamFunction<double> f;
	bcc_odd<quintic_box<double,double>, double, double> lat(1./double(2.*100));

//	viennacl::ocl::set_context_device_type(0, viennacl::ocl::cpu_tag());
//	printf("OpenCL Info:\n%s \n", viennacl::ocl::current_device().info().c_str());

	// Make a pointset
	std::vector<vector3<double>> pts(N_VAR_SAMPLES);
	std::vector<double> smpls(N_VAR_SAMPLES);

	#pragma omp parallel for
	for(int i = 0; i < N_VAR_SAMPLES; i++){
		vector3<double>  p(	((double)rand()/(double)RAND_MAX),
							((double)rand()/(double)RAND_MAX),
							((double)rand()/(double)RAND_MAX));
		pts[i] = p;
		smpls[i] = f.evaluate(p);
	}
	std::cout << sisl::utility::_bcc_optimization_kernel_dbl << std :: endl;

	utility::splatSamplesOntoLattice<lType, gType, double, double>(
				pts,
				smpls,
				&lat,
				vector3<double>(0,0,0));

	utility::fitPoints<lType, gType, double, double> (
				pts,
				smpls,
				gType::getBeppoLevi2Norm(),
				&lat,
				vector3<double>(0,0,0)
			);

	mc.marchLattice<lType, double, double>(
		&lat, 
		NULL, 
		NULL, 
		NULL, 
		0.5, 
		0.01, 
		vector3<double>(0,0,0), 
		vector3<double>(1,1,1));
	
	REQUIRE(mc.writeSurface("scattered_bcc_test.ply"));
}
