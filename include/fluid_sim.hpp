#include <sisl/sisl.hpp>
#include <sisl/lattice/cartesian_odd.hpp>
#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/basis/tp3cubic.hpp>
#include <sisl/basis/quintic.hpp>

using namespace sisl;

template<class I, class O, class L>
class fluid_sim {
public:
	fluid_sim(I h, I stepSize){

	}

	void setForceFunction(){

	}

	void step(){
		advect();
		diffuse();
		addExternalForces();
		project();
	}

	O eval(I x, I y, I z){

	}
private:
	void advect();
	void diffuse();
	void addExternalForces();
	void project();
};

template<class I, class O>
class fluid_sim<I, O, bcc_odd<quintic_box<I,O>, I, O> > {
public:
	fluid_sim(I h, I stepSize){

	}

	void setForceFunction(){

	}

	void step(){
		advect();
		diffuse();
		addExternalForces();
		project();
	}

	O eval(I x, I y, I z){

	}
private:
	void advect() {

	}
	void diffuse();
	void addExternalForces();
	void project();
	bcc_odd<quintic_box<I,O>, I, O> u;
	bcc_odd<quintic_box<I,O>, I, O> v;
	bcc_odd<quintic_box<I,O>, I, O> w;
};

template<class I, class O>
class fluid_sim<I, O, cartesian_odd<tp3cubic<I,O>, I, O> > {
public:
	fluid_sim(I h, I stepSize){
	}

	void setForceFunction(){

	}

	void step(){
		advect();
		diffuse();
		addExternalForces();
		project();
	}

	O eval(I x, I y, I z){

	}
private:
	void advect() {
		// U sim
		u.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			vector3<double> p = u.getSitePosition(i,j,k);
			return p.i - timestep * u.f(p);
		});

		v.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			vector3<double> p = v.getSitePosition(i,j,k);
			return p.j - timestep * v.f(p);
		});

		w.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			vector3<double> p = w.getSitePosition(i,j,k);
			return p.k - timestep * w.f(p);
		});

	}
	void diffuse(){}
	void addExternalForces(){}
	void project(){
		auto dx = tp3cubic<I,O>::dxFilter0bd(&u);
		auto dy = tp3cubic<I,O>::dxFilter0bd(&v);
		auto dz = tp3cubic<I,O>::dxFilter0bd(&w);
		auto pAccess  = 
			[&](const int &x, const int &y, const int &z) ->  O {
				return (O) p.GV(x,y,z); 
			};
		
		cartesian_odd<tp3cubic<I,O>, I, O> *ur = &u;
		cartesian_odd<tp3cubic<I,O>, I, O> *vr = &v;
		cartesian_odd<tp3cubic<I,O>, I, O> *wr = &w;

		p.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			return u.getSitePosition(i,j,k);
		});

		// Take divergenece of velocity field
		p.forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double sx = 1./(x->getScale()*12.);

			double dx = (
			 -1.*ur->GV(i-2, j, k) + 
			  8.*ur->GV(i-1, j, k) +
			 -8.*ur->GV(i+1, j, k) + 
				 ur->GV(i+2, j, k)) * sx;
		
			double dy = (
			 (-1.*vr->GV(i, j-2, k)) + 
			 ( 8.*vr->GV(i, j-1, k)) + 
			 (-8.*vr->GV(i, j+1, k)) + 
			 (	vr->GV(i, j+2, k))) * sx;

			double dz = (
			 -1.*wr->GV(i, j, k-2) + 
			  8.*wr->GV(i, j, k-1) + 
			 -8.*wr->GV(i, j, k+1) + 
				 wr->GV(i, j, k+2)) * sx; 

			return (dx+dy+dz);
		});


		// Solve the poisson problem
		p.frequencyFilter(
			combineFilters<double,double>(
				tp3cubic<double,double>::interpolationFilter,
				tp3cubic<double,double>::poissonFilter
			)
		);


		// Subtract the pressure gradient from the vector feild
		u.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			return u.GV(i,j,k) - dx(i,j,k,pAccess);
		});
		v.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			return v.GV(i,j,k) - dy(i,j,k,pAccess);
		});
		w.spatialFilter([&](const int &i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s) {
			return w.GV(i,j,k) - dz(i,j,k,pAccess);
		});

		u.frequencyFilter(tp3cubic<double,double>::interpolationFilter);
		v.frequencyFilter(tp3cubic<double,double>::interpolationFilter);
		w.frequencyFilter(tp3cubic<double,double>::interpolationFilter);
	}

	I timestep;
	cartesian_odd<tp3cubic<I,O>, I, O> u;
	cartesian_odd<tp3cubic<I,O>, I, O> v;
	cartesian_odd<tp3cubic<I,O>, I, O> w;
	cartesian_odd<tp3cubic<I,O>, I, O> p;
};