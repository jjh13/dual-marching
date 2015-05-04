#define CATCH_CONFIG_MAIN 
#include <catch/catch.hpp>

#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/basis/quintic.hpp>

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


TEST_CASE("Spatial filter test", "sfilt"){
	typedef bcc_odd<quintic_box<double,double>, double, double>  LATType;

	LATType *x = new LATType(1./double(2.*100));
	LATType *y = new LATType(1./double(2.*100));
	LATType *z = new LATType(1./double(2.*100));
	LATType *l = new LATType(1./double(2.*100));
	
	std::function<double(const int &, const int &, const int &, std::function<double(const int &, const int &, const int & )>)> 
		divFunction;


	l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double iscl = -1./(x->getScale()*12.); 

			double dx = 
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

			return iscl*(dx+dy+dz);
		} 
	);
}

TEST_CASE("Array Write Test","arraywrite"){
	using namespace sisl;
	TestFunction<double> f(6.);
	int idx = 0;
	double error[1000];

	for(int res = 25; res <= 200; res = res * 2, idx++){
		bcc_odd<quintic_box<double,double>, double, double> lat(1./double(2.*res));
		error[idx] = 0;

		printf("Sampling lattice of size (%d*%d*%d)...\n", res, res, res);

		lat.forEachLatticeSite([&](const int &i, const int &j, const int &k) {
			vector3<double> p = lat.getSitePosition(i,j,k);
			return f.evaluate(p);
		});

		printf("Prefilering...\n");	
		lat.frequencyFilter(quintic_box<double,double>::interpolationFilter);

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

TEST_CASE("Write Test","arraywrite"){
	using namespace sisl;
	TestFunction<double> f(6.);
	int idx = 0;
	double error[1000];

	for(int res = 25; res <= 200; res = res * 2, idx++){
		bcc_odd<quintic_box<double,double>, double, double> lat(1./double(2.*res));
		error[idx] = 0;

		printf("Sampling lattice of size (%d*%d*%d)...\n", res, res, res);

		lat.forEachLatticeSite([&](const int &i, const int &j, const int &k) {
			vector3<double> p = lat.getSitePosition(i,j,k);
			return f.laplacian(p);
		});

		printf("Prefilering...\n");	
		lat.frequencyFilter(
			combineFilters<double,double>(
				quintic_box<double,double>::interpolationFilter,
				quintic_box<double,double>::poissonFilter
			)
		);

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
