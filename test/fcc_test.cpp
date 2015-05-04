#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>

#include <sisl/lattice/fcc_odd.hpp>
#include <sisl/basis/delta_basis.hpp>

#include <iostream>

#define N_MC_SAMPLES  100000
#define N_VAR_SAMPLES 300000

using namespace sisl;

template<class I>
class TestFunction {
	I mu;
public:
	TestFunction(I mew) : mu(mew){}
	//! Evaluate a function at a point
	I evaluate(const I &x, const I &y, const I &z) {
		return(sin(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)));
	}
	I evaluate(const vector3<I> &p) {
		return this->evaluate(p.i, p.j, p.k);
	}
	
	I div(const vector3<I> &p) { return this->div(p.i, p.j, p.k); }
	I div(const I &x, const I &y, const I &z) {
		return 0.2e1 * cos(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * 0.3141592654e1 * 0.3141592654e1 * mu * cos(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1) + 0.2e1 * cos(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * 0.3141592654e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * cos(y * 0.3141592654e1) * sin(z * 0.3141592654e1) + 0.2e1 * cos(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * 0.3141592654e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * cos(z * 0.3141592654e1);
	}
	
	I laplacian(const vector3<I> &p){return this->laplacian(p.i, p.j, p.k);}
	I laplacian(const I &x, const I &y, const I &z) {
		return -0.4e1 * sin(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * pow(0.3141592654e1, 0.4e1) * mu * mu * pow(cos(x * 0.3141592654e1), 0.2e1) * pow(sin(y * 0.3141592654e1), 0.2e1) * pow(sin(z * 0.3141592654e1), 0.2e1) - 0.6e1 * cos(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * pow(0.3141592654e1, 0.3e1) * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1) - 0.4e1 * sin(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * pow(0.3141592654e1, 0.4e1) * mu * mu * pow(sin(x * 0.3141592654e1), 0.2e1) * pow(cos(y * 0.3141592654e1), 0.2e1) * pow(sin(z * 0.3141592654e1), 0.2e1) - 0.4e1 * sin(0.2e1 * 0.3141592654e1 * mu * sin(x * 0.3141592654e1) * sin(y * 0.3141592654e1) * sin(z * 0.3141592654e1)) * pow(0.3141592654e1, 0.4e1) * mu * mu * pow(sin(x * 0.3141592654e1), 0.2e1) * pow(sin(y * 0.3141592654e1), 0.2e1) * pow(cos(z * 0.3141592654e1), 0.2e1);
	}
};


TEST_CASE("Array Write Test","arraywrite"){
	using namespace sisl;
	TestFunction<double> f(6.);
	int idx = 0;
	double error[1000];

	for(int res = 32; res <= 256; res = res * 2, idx++){
		fcc_odd<delta_basis<double,double>, double, double> lat(1./double(res-1));
		error[idx] = 0;

		printf("Sampling lattice of size (%d*%d*%d)...\n", res, res, res);

		#pragma omp parallel for
		for(int i = 1; i < res - 1; i++)
			for(int j = 1; j < res - 1; j++)
				for(int k = 1; k < res - 1; k++){
					vector3<double> p = lat.getSitePosition(i,j,k);
					lat.SV(i,j,k, f.evaluate(p));
				}

		printf("Prefilering...\n");	
		//lat.frequencyFilter(tp3cubic<double,double>::interpolationFilter);

		printf("Reonconstruction test...\n");
		#pragma omp parallel for
		for(int i = 0; i < N_MC_SAMPLES; i++){
       		vector3<double>  p(	((double)rand()/(double)RAND_MAX),
       					((double)rand()/(double)RAND_MAX),
       					((double)rand()/(double)RAND_MAX));
       		double t = powf(lat.f(p) - f.evaluate(p), 2.);
       		#pragma omp atomic
        	error[idx] += t;
		}
		error[idx] /= double(N_MC_SAMPLES);
		error[idx] = sqrt(error[idx]);
		printf("msl2e = %e\n", error[idx]); 
		if(idx > 0)
			printf("Log drop from last step: %f\n", (log(error[idx-1]/error[idx]))/log(2.) );
	}
}