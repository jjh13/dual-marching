#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _TP_BSPLINE3_H_
#define _TP_BSPLINE3_H_

namespace sisl{
template <class O, class I>
class tp3cubic: public basis_function <O,I> {
private:
	static inline O bspline3(const I &t){
		if(t <= -2  || t > 2)
			return 0;
		if(t <= -1 && t > -2)
			return (4./3.) + t*(2. + t*(1. + (1./6.)*t));
		if(t <= 0 && t > -1)
			return (2./3.) + t*t*(-1. - 0.5*t);
		if( t <= 1 && t > 0)
			return (2./3.) + t*t*(-1. + 0.5*t);
		return (4./3.) + t*(-2. + t*(1. - (1./6.)*t));
	}

public:
	static std::string getBasisName(){
		return std::string("tp3cubic");
	}
	static const O M(const I &x, const I &y, const I &z) { return (O)bspline3(x)*bspline3(y)*bspline3(z); };
	static const O M(const vector3<I> &p) { return (O)bspline3(p.i)*bspline3(p.j)*bspline3(p.k); };

	// This function should return the intersection of the closure of the support of
	// the generator and the lattice.
	static std::vector<std::tuple<int,int,int>> getSupport() {
		using namespace std;
		std::vector<std::tuple<int,int,int>> support;
		for(int i = -2; i <= 2; i++)
			for(int j = -2; j <= 2; j++)
				for(int k = -2; k <= 2; k++){
					support.push_back(make_tuple(i,j,k));
				}
		return support;
	};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int>> getEffectiveSupport(const vector3<I> &p) {
		using namespace std;
		std::vector<std::tuple<int,int,int>> support;
		for(int i = -2; i <= 2; i++)
			for(int j = -2; j <= 2; j++)
				for(int k = -2; k <= 2; k++){
					support.push_back(make_tuple(i,j,k));
				}
		return support;
	}

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi2Norm(){
		using namespace std;
		return {make_tuple(-3, -3, -3, 1.0235575711766190e-07), make_tuple(-3, -3, -2, 6.2043335852859660e-06), make_tuple(-3, -3, -1, 4.9217372134038800e-05), make_tuple(-3, -3, 0, 9.3978332073570170e-05), make_tuple(-3, -3, 1, 4.9217372134038800e-05), make_tuple(-3, -3, 2, 6.2043335852859660e-06), make_tuple(-3, -3, 3, 1.0235575711766190e-07), make_tuple(-3, -2, -3, 6.2043335852859660e-06), make_tuple(-3, -2, -2, 2.6908541194255480e-04), make_tuple(-3, -2, -1, 1.7778722600151170e-03), make_tuple(-3, -2, 0, 3.1952632905013860e-03), make_tuple(-3, -2, 1, 1.7778722600151170e-03), make_tuple(-3, -2, 2, 2.6908541194255480e-04), make_tuple(-3, -2, 3, 6.2043335852859660e-06), make_tuple(-3, -1, -3, 4.9217372134038800e-05), make_tuple(-3, -1, -2, 1.7778722600151170e-03), make_tuple(-3, -1, -1, 1.0157194822373390e-02), make_tuple(-3, -1, 0, 1.7202066011589820e-02), make_tuple(-3, -1, 1, 1.0157194822373390e-02), make_tuple(-3, -1, 2, 1.7778722600151170e-03), make_tuple(-3, -1, 3, 4.9217372134038800e-05), make_tuple(-3, 0, -3, 9.3978332073570170e-05), make_tuple(-3, 0, -2, 3.1952632905013860e-03), make_tuple(-3, 0, -1, 1.7202066011589820e-02), make_tuple(-3, 0, 0, 2.8329554043839760e-02), make_tuple(-3, 0, 1, 1.7202066011589820e-02), make_tuple(-3, 0, 2, 3.1952632905013860e-03), make_tuple(-3, 0, 3, 9.3978332073570170e-05), make_tuple(-3, 1, -3, 4.9217372134038800e-05), make_tuple(-3, 1, -2, 1.7778722600151170e-03), make_tuple(-3, 1, -1, 1.0157194822373390e-02), make_tuple(-3, 1, 0, 1.7202066011589820e-02), make_tuple(-3, 1, 1, 1.0157194822373390e-02), make_tuple(-3, 1, 2, 1.7778722600151170e-03), make_tuple(-3, 1, 3, 4.9217372134038800e-05), make_tuple(-3, 2, -3, 6.2043335852859660e-06), make_tuple(-3, 2, -2, 2.6908541194255480e-04), make_tuple(-3, 2, -1, 1.7778722600151170e-03), make_tuple(-3, 2, 0, 3.1952632905013860e-03), make_tuple(-3, 2, 1, 1.7778722600151170e-03), make_tuple(-3, 2, 2, 2.6908541194255480e-04), make_tuple(-3, 2, 3, 6.2043335852859660e-06), make_tuple(-3, 3, -3, 1.0235575711766190e-07), make_tuple(-3, 3, -2, 6.2043335852859660e-06), make_tuple(-3, 3, -1, 4.9217372134038800e-05), make_tuple(-3, 3, 0, 9.3978332073570170e-05), make_tuple(-3, 3, 1, 4.9217372134038800e-05), make_tuple(-3, 3, 2, 6.2043335852859660e-06), make_tuple(-3, 3, 3, 1.0235575711766190e-07), make_tuple(-2, -3, -3, 6.2043335852859660e-06), make_tuple(-2, -3, -2, 2.6908541194255480e-04), make_tuple(-2, -3, -1, 1.7778722600151170e-03), make_tuple(-2, -3, 0, 3.1952632905013860e-03), make_tuple(-2, -3, 1, 1.7778722600151170e-03), make_tuple(-2, -3, 2, 2.6908541194255480e-04), make_tuple(-2, -3, 3, 6.2043335852859660e-06), make_tuple(-2, -2, -3, 2.6908541194255480e-04), make_tuple(-2, -2, -2, 5.7142857142857140e-03), make_tuple(-2, -2, -1, 2.0435374149659860e-02), make_tuple(-2, -2, 0, 2.7162509448223730e-02), make_tuple(-2, -2, 1, 2.0435374149659860e-02), make_tuple(-2, -2, 2, 5.7142857142857140e-03), make_tuple(-2, -2, 3, 2.6908541194255480e-04), make_tuple(-2, -1, -3, 1.7778722600151170e-03), make_tuple(-2, -1, -2, 2.0435374149659860e-02), make_tuple(-2, -1, -1, 7.4957482993197280e-03), make_tuple(-2, -1, 0, -4.5132275132275130e-02), make_tuple(-2, -1, 1, 7.4957482993197280e-03), make_tuple(-2, -1, 2, 2.0435374149659860e-02), make_tuple(-2, -1, 3, 1.7778722600151170e-03), make_tuple(-2, 0, -3, 3.1952632905013860e-03), make_tuple(-2, 0, -2, 2.7162509448223730e-02), make_tuple(-2, 0, -1, -4.5132275132275130e-02), make_tuple(-2, 0, 0, -1.7362559838750310e-01), make_tuple(-2, 0, 1, -4.5132275132275130e-02), make_tuple(-2, 0, 2, 2.7162509448223730e-02), make_tuple(-2, 0, 3, 3.1952632905013860e-03), make_tuple(-2, 1, -3, 1.7778722600151170e-03), make_tuple(-2, 1, -2, 2.0435374149659860e-02), make_tuple(-2, 1, -1, 7.4957482993197280e-03), make_tuple(-2, 1, 0, -4.5132275132275130e-02), make_tuple(-2, 1, 1, 7.4957482993197280e-03), make_tuple(-2, 1, 2, 2.0435374149659860e-02), make_tuple(-2, 1, 3, 1.7778722600151170e-03), make_tuple(-2, 2, -3, 2.6908541194255480e-04), make_tuple(-2, 2, -2, 5.7142857142857140e-03), make_tuple(-2, 2, -1, 2.0435374149659860e-02), make_tuple(-2, 2, 0, 2.7162509448223730e-02), make_tuple(-2, 2, 1, 2.0435374149659860e-02), make_tuple(-2, 2, 2, 5.7142857142857140e-03), make_tuple(-2, 2, 3, 2.6908541194255480e-04), make_tuple(-2, 3, -3, 6.2043335852859660e-06), make_tuple(-2, 3, -2, 2.6908541194255480e-04), make_tuple(-2, 3, -1, 1.7778722600151170e-03), make_tuple(-2, 3, 0, 3.1952632905013860e-03), make_tuple(-2, 3, 1, 1.7778722600151170e-03), make_tuple(-2, 3, 2, 2.6908541194255480e-04), make_tuple(-2, 3, 3, 6.2043335852859660e-06), make_tuple(-1, -3, -3, 4.9217372134038800e-05), make_tuple(-1, -3, -2, 1.7778722600151170e-03), make_tuple(-1, -3, -1, 1.0157194822373390e-02), make_tuple(-1, -3, 0, 1.7202066011589820e-02), make_tuple(-1, -3, 1, 1.0157194822373390e-02), make_tuple(-1, -3, 2, 1.7778722600151170e-03), make_tuple(-1, -3, 3, 4.9217372134038800e-05), make_tuple(-1, -2, -3, 1.7778722600151170e-03), make_tuple(-1, -2, -2, 2.0435374149659860e-02), make_tuple(-1, -2, -1, 7.4957482993197280e-03), make_tuple(-1, -2, 0, -4.5132275132275130e-02), make_tuple(-1, -2, 1, 7.4957482993197280e-03), make_tuple(-1, -2, 2, 2.0435374149659860e-02), make_tuple(-1, -2, 3, 1.7778722600151170e-03), make_tuple(-1, -1, -3, 1.0157194822373390e-02), make_tuple(-1, -1, -2, 7.4957482993197280e-03), make_tuple(-1, -1, -1, -2.2913584183673470e-01), make_tuple(-1, -1, 0, -2.5471277399848830e-01), make_tuple(-1, -1, 1, -2.2913584183673470e-01), make_tuple(-1, -1, 2, 7.4957482993197280e-03), make_tuple(-1, -1, 3, 1.0157194822373390e-02), make_tuple(-1, 0, -3, 1.7202066011589820e-02), make_tuple(-1, 0, -2, -4.5132275132275130e-02), make_tuple(-1, 0, -1, -2.5471277399848830e-01), make_tuple(-1, 0, 0, 3.0973041068279160e-01), make_tuple(-1, 0, 1, -2.5471277399848830e-01), make_tuple(-1, 0, 2, -4.5132275132275130e-02), make_tuple(-1, 0, 3, 1.7202066011589820e-02), make_tuple(-1, 1, -3, 1.0157194822373390e-02), make_tuple(-1, 1, -2, 7.4957482993197280e-03), make_tuple(-1, 1, -1, -2.2913584183673470e-01), make_tuple(-1, 1, 0, -2.5471277399848830e-01), make_tuple(-1, 1, 1, -2.2913584183673470e-01), make_tuple(-1, 1, 2, 7.4957482993197280e-03), make_tuple(-1, 1, 3, 1.0157194822373390e-02), make_tuple(-1, 2, -3, 1.7778722600151170e-03), make_tuple(-1, 2, -2, 2.0435374149659860e-02), make_tuple(-1, 2, -1, 7.4957482993197280e-03), make_tuple(-1, 2, 0, -4.5132275132275130e-02), make_tuple(-1, 2, 1, 7.4957482993197280e-03), make_tuple(-1, 2, 2, 2.0435374149659860e-02), make_tuple(-1, 2, 3, 1.7778722600151170e-03), make_tuple(-1, 3, -3, 4.9217372134038800e-05), make_tuple(-1, 3, -2, 1.7778722600151170e-03), make_tuple(-1, 3, -1, 1.0157194822373390e-02), make_tuple(-1, 3, 0, 1.7202066011589820e-02), make_tuple(-1, 3, 1, 1.0157194822373390e-02), make_tuple(-1, 3, 2, 1.7778722600151170e-03), make_tuple(-1, 3, 3, 4.9217372134038800e-05), make_tuple(0, -3, -3, 9.3978332073570170e-05), make_tuple(0, -3, -2, 3.1952632905013860e-03), make_tuple(0, -3, -1, 1.7202066011589820e-02), make_tuple(0, -3, 0, 2.8329554043839760e-02), make_tuple(0, -3, 1, 1.7202066011589820e-02), make_tuple(0, -3, 2, 3.1952632905013860e-03), make_tuple(0, -3, 3, 9.3978332073570170e-05), make_tuple(0, -2, -3, 3.1952632905013860e-03), make_tuple(0, -2, -2, 2.7162509448223730e-02), make_tuple(0, -2, -1, -4.5132275132275130e-02), make_tuple(0, -2, 0, -1.7362559838750310e-01), make_tuple(0, -2, 1, -4.5132275132275130e-02), make_tuple(0, -2, 2, 2.7162509448223730e-02), make_tuple(0, -2, 3, 3.1952632905013860e-03), make_tuple(0, -1, -3, 1.7202066011589820e-02), make_tuple(0, -1, -2, -4.5132275132275130e-02), make_tuple(0, -1, -1, -2.5471277399848830e-01), make_tuple(0, -1, 0, 3.0973041068279160e-01), make_tuple(0, -1, 1, -2.5471277399848830e-01), make_tuple(0, -1, 2, -4.5132275132275130e-02), make_tuple(0, -1, 3, 1.7202066011589820e-02), make_tuple(0, 0, -3, 2.8329554043839760e-02), make_tuple(0, 0, -2, -1.7362559838750310e-01), make_tuple(0, 0, -1, 3.0973041068279160e-01), make_tuple(0, 0, 0, 3.1166339128243890e+00), make_tuple(0, 0, 1, 3.0973041068279160e-01), make_tuple(0, 0, 2, -1.7362559838750310e-01), make_tuple(0, 0, 3, 2.8329554043839760e-02), make_tuple(0, 1, -3, 1.7202066011589820e-02), make_tuple(0, 1, -2, -4.5132275132275130e-02), make_tuple(0, 1, -1, -2.5471277399848830e-01), make_tuple(0, 1, 0, 3.0973041068279160e-01), make_tuple(0, 1, 1, -2.5471277399848830e-01), make_tuple(0, 1, 2, -4.5132275132275130e-02), make_tuple(0, 1, 3, 1.7202066011589820e-02), make_tuple(0, 2, -3, 3.1952632905013860e-03), make_tuple(0, 2, -2, 2.7162509448223730e-02), make_tuple(0, 2, -1, -4.5132275132275130e-02), make_tuple(0, 2, 0, -1.7362559838750310e-01), make_tuple(0, 2, 1, -4.5132275132275130e-02), make_tuple(0, 2, 2, 2.7162509448223730e-02), make_tuple(0, 2, 3, 3.1952632905013860e-03), make_tuple(0, 3, -3, 9.3978332073570170e-05), make_tuple(0, 3, -2, 3.1952632905013860e-03), make_tuple(0, 3, -1, 1.7202066011589820e-02), make_tuple(0, 3, 0, 2.8329554043839760e-02), make_tuple(0, 3, 1, 1.7202066011589820e-02), make_tuple(0, 3, 2, 3.1952632905013860e-03), make_tuple(0, 3, 3, 9.3978332073570170e-05), make_tuple(1, -3, -3, 4.9217372134038800e-05), make_tuple(1, -3, -2, 1.7778722600151170e-03), make_tuple(1, -3, -1, 1.0157194822373390e-02), make_tuple(1, -3, 0, 1.7202066011589820e-02), make_tuple(1, -3, 1, 1.0157194822373390e-02), make_tuple(1, -3, 2, 1.7778722600151170e-03), make_tuple(1, -3, 3, 4.9217372134038800e-05), make_tuple(1, -2, -3, 1.7778722600151170e-03), make_tuple(1, -2, -2, 2.0435374149659860e-02), make_tuple(1, -2, -1, 7.4957482993197280e-03), make_tuple(1, -2, 0, -4.5132275132275130e-02), make_tuple(1, -2, 1, 7.4957482993197280e-03), make_tuple(1, -2, 2, 2.0435374149659860e-02), make_tuple(1, -2, 3, 1.7778722600151170e-03), make_tuple(1, -1, -3, 1.0157194822373390e-02), make_tuple(1, -1, -2, 7.4957482993197280e-03), make_tuple(1, -1, -1, -2.2913584183673470e-01), make_tuple(1, -1, 0, -2.5471277399848830e-01), make_tuple(1, -1, 1, -2.2913584183673470e-01), make_tuple(1, -1, 2, 7.4957482993197280e-03), make_tuple(1, -1, 3, 1.0157194822373390e-02), make_tuple(1, 0, -3, 1.7202066011589820e-02), make_tuple(1, 0, -2, -4.5132275132275130e-02), make_tuple(1, 0, -1, -2.5471277399848830e-01), make_tuple(1, 0, 0, 3.0973041068279160e-01), make_tuple(1, 0, 1, -2.5471277399848830e-01), make_tuple(1, 0, 2, -4.5132275132275130e-02), make_tuple(1, 0, 3, 1.7202066011589820e-02), make_tuple(1, 1, -3, 1.0157194822373390e-02), make_tuple(1, 1, -2, 7.4957482993197280e-03), make_tuple(1, 1, -1, -2.2913584183673470e-01), make_tuple(1, 1, 0, -2.5471277399848830e-01), make_tuple(1, 1, 1, -2.2913584183673470e-01), make_tuple(1, 1, 2, 7.4957482993197280e-03), make_tuple(1, 1, 3, 1.0157194822373390e-02), make_tuple(1, 2, -3, 1.7778722600151170e-03), make_tuple(1, 2, -2, 2.0435374149659860e-02), make_tuple(1, 2, -1, 7.4957482993197280e-03), make_tuple(1, 2, 0, -4.5132275132275130e-02), make_tuple(1, 2, 1, 7.4957482993197280e-03), make_tuple(1, 2, 2, 2.0435374149659860e-02), make_tuple(1, 2, 3, 1.7778722600151170e-03), make_tuple(1, 3, -3, 4.9217372134038800e-05), make_tuple(1, 3, -2, 1.7778722600151170e-03), make_tuple(1, 3, -1, 1.0157194822373390e-02), make_tuple(1, 3, 0, 1.7202066011589820e-02), make_tuple(1, 3, 1, 1.0157194822373390e-02), make_tuple(1, 3, 2, 1.7778722600151170e-03), make_tuple(1, 3, 3, 4.9217372134038800e-05), make_tuple(2, -3, -3, 6.2043335852859660e-06), make_tuple(2, -3, -2, 2.6908541194255480e-04), make_tuple(2, -3, -1, 1.7778722600151170e-03), make_tuple(2, -3, 0, 3.1952632905013860e-03), make_tuple(2, -3, 1, 1.7778722600151170e-03), make_tuple(2, -3, 2, 2.6908541194255480e-04), make_tuple(2, -3, 3, 6.2043335852859660e-06), make_tuple(2, -2, -3, 2.6908541194255480e-04), make_tuple(2, -2, -2, 5.7142857142857140e-03), make_tuple(2, -2, -1, 2.0435374149659860e-02), make_tuple(2, -2, 0, 2.7162509448223730e-02), make_tuple(2, -2, 1, 2.0435374149659860e-02), make_tuple(2, -2, 2, 5.7142857142857140e-03), make_tuple(2, -2, 3, 2.6908541194255480e-04), make_tuple(2, -1, -3, 1.7778722600151170e-03), make_tuple(2, -1, -2, 2.0435374149659860e-02), make_tuple(2, -1, -1, 7.4957482993197280e-03), make_tuple(2, -1, 0, -4.5132275132275130e-02), make_tuple(2, -1, 1, 7.4957482993197280e-03), make_tuple(2, -1, 2, 2.0435374149659860e-02), make_tuple(2, -1, 3, 1.7778722600151170e-03), make_tuple(2, 0, -3, 3.1952632905013860e-03), make_tuple(2, 0, -2, 2.7162509448223730e-02), make_tuple(2, 0, -1, -4.5132275132275130e-02), make_tuple(2, 0, 0, -1.7362559838750310e-01), make_tuple(2, 0, 1, -4.5132275132275130e-02), make_tuple(2, 0, 2, 2.7162509448223730e-02), make_tuple(2, 0, 3, 3.1952632905013860e-03), make_tuple(2, 1, -3, 1.7778722600151170e-03), make_tuple(2, 1, -2, 2.0435374149659860e-02), make_tuple(2, 1, -1, 7.4957482993197280e-03), make_tuple(2, 1, 0, -4.5132275132275130e-02), make_tuple(2, 1, 1, 7.4957482993197280e-03), make_tuple(2, 1, 2, 2.0435374149659860e-02), make_tuple(2, 1, 3, 1.7778722600151170e-03), make_tuple(2, 2, -3, 2.6908541194255480e-04), make_tuple(2, 2, -2, 5.7142857142857140e-03), make_tuple(2, 2, -1, 2.0435374149659860e-02), make_tuple(2, 2, 0, 2.7162509448223730e-02), make_tuple(2, 2, 1, 2.0435374149659860e-02), make_tuple(2, 2, 2, 5.7142857142857140e-03), make_tuple(2, 2, 3, 2.6908541194255480e-04), make_tuple(2, 3, -3, 6.2043335852859660e-06), make_tuple(2, 3, -2, 2.6908541194255480e-04), make_tuple(2, 3, -1, 1.7778722600151170e-03), make_tuple(2, 3, 0, 3.1952632905013860e-03), make_tuple(2, 3, 1, 1.7778722600151170e-03), make_tuple(2, 3, 2, 2.6908541194255480e-04), make_tuple(2, 3, 3, 6.2043335852859660e-06), make_tuple(3, -3, -3, 1.0235575711766190e-07), make_tuple(3, -3, -2, 6.2043335852859660e-06), make_tuple(3, -3, -1, 4.9217372134038800e-05), make_tuple(3, -3, 0, 9.3978332073570170e-05), make_tuple(3, -3, 1, 4.9217372134038800e-05), make_tuple(3, -3, 2, 6.2043335852859660e-06), make_tuple(3, -3, 3, 1.0235575711766190e-07), make_tuple(3, -2, -3, 6.2043335852859660e-06), make_tuple(3, -2, -2, 2.6908541194255480e-04), make_tuple(3, -2, -1, 1.7778722600151170e-03), make_tuple(3, -2, 0, 3.1952632905013860e-03), make_tuple(3, -2, 1, 1.7778722600151170e-03), make_tuple(3, -2, 2, 2.6908541194255480e-04), make_tuple(3, -2, 3, 6.2043335852859660e-06), make_tuple(3, -1, -3, 4.9217372134038800e-05), make_tuple(3, -1, -2, 1.7778722600151170e-03), make_tuple(3, -1, -1, 1.0157194822373390e-02), make_tuple(3, -1, 0, 1.7202066011589820e-02), make_tuple(3, -1, 1, 1.0157194822373390e-02), make_tuple(3, -1, 2, 1.7778722600151170e-03), make_tuple(3, -1, 3, 4.9217372134038800e-05), make_tuple(3, 0, -3, 9.3978332073570170e-05), make_tuple(3, 0, -2, 3.1952632905013860e-03), make_tuple(3, 0, -1, 1.7202066011589820e-02), make_tuple(3, 0, 0, 2.8329554043839760e-02), make_tuple(3, 0, 1, 1.7202066011589820e-02), make_tuple(3, 0, 2, 3.1952632905013860e-03), make_tuple(3, 0, 3, 9.3978332073570170e-05), make_tuple(3, 1, -3, 4.9217372134038800e-05), make_tuple(3, 1, -2, 1.7778722600151170e-03), make_tuple(3, 1, -1, 1.0157194822373390e-02), make_tuple(3, 1, 0, 1.7202066011589820e-02), make_tuple(3, 1, 1, 1.0157194822373390e-02), make_tuple(3, 1, 2, 1.7778722600151170e-03), make_tuple(3, 1, 3, 4.9217372134038800e-05), make_tuple(3, 2, -3, 6.2043335852859660e-06), make_tuple(3, 2, -2, 2.6908541194255480e-04), make_tuple(3, 2, -1, 1.7778722600151170e-03), make_tuple(3, 2, 0, 3.1952632905013860e-03), make_tuple(3, 2, 1, 1.7778722600151170e-03), make_tuple(3, 2, 2, 2.6908541194255480e-04), make_tuple(3, 2, 3, 6.2043335852859660e-06), make_tuple(3, 3, -3, 1.0235575711766190e-07), make_tuple(3, 3, -2, 6.2043335852859660e-06), make_tuple(3, 3, -1, 4.9217372134038800e-05), make_tuple(3, 3, 0, 9.3978332073570170e-05), make_tuple(3, 3, 1, 4.9217372134038800e-05), make_tuple(3, 3, 2, 6.2043335852859660e-06), make_tuple(3, 3, 3, 1.0235575711766190e-07)};
	};
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> autoCorrelation(){
		using namespace std;
		return {make_tuple(-3, -3, -3, 7.8110315260730980e-12), make_tuple(-3, -3, -2, 9.3732378312877180e-10), make_tuple(-3, -3, -1, 9.3029385475530600e-09), make_tuple(-3, -3, 0, 1.8871452166992610e-08), make_tuple(-3, -3, 1, 9.3029385475530600e-09), make_tuple(-3, -3, 2, 9.3732378312877180e-10), make_tuple(-3, -3, 3, 7.8110315260730980e-12), make_tuple(-3, -2, -3, 9.3732378312877180e-10), make_tuple(-3, -2, -2, 1.1247885397545260e-07), make_tuple(-3, -2, -1, 1.1163526257063670e-06), make_tuple(-3, -2, 0, 2.2645742600391130e-06), make_tuple(-3, -2, 1, 1.1163526257063670e-06), make_tuple(-3, -2, 2, 1.1247885397545260e-07), make_tuple(-3, -2, 3, 9.3732378312877180e-10), make_tuple(-3, -1, -3, 9.3029385475530600e-09), make_tuple(-3, -1, -2, 1.1163526257063670e-06), make_tuple(-3, -1, -1, 1.1079799810135690e-05), make_tuple(-3, -1, 0, 2.2475899530888190e-05), make_tuple(-3, -1, 1, 1.1079799810135690e-05), make_tuple(-3, -1, 2, 1.1163526257063670e-06), make_tuple(-3, -1, 3, 9.3029385475530600e-09), make_tuple(-3, 0, -3, 1.8871452166992610e-08), make_tuple(-3, 0, -2, 2.2645742600391130e-06), make_tuple(-3, 0, -1, 2.2475899530888190e-05), make_tuple(-3, 0, 0, 4.5593428435454130e-05), make_tuple(-3, 0, 1, 2.2475899530888190e-05), make_tuple(-3, 0, 2, 2.2645742600391130e-06), make_tuple(-3, 0, 3, 1.8871452166992610e-08), make_tuple(-3, 1, -3, 9.3029385475530600e-09), make_tuple(-3, 1, -2, 1.1163526257063670e-06), make_tuple(-3, 1, -1, 1.1079799810135690e-05), make_tuple(-3, 1, 0, 2.2475899530888190e-05), make_tuple(-3, 1, 1, 1.1079799810135690e-05), make_tuple(-3, 1, 2, 1.1163526257063670e-06), make_tuple(-3, 1, 3, 9.3029385475530600e-09), make_tuple(-3, 2, -3, 9.3732378312877180e-10), make_tuple(-3, 2, -2, 1.1247885397545260e-07), make_tuple(-3, 2, -1, 1.1163526257063670e-06), make_tuple(-3, 2, 0, 2.2645742600391130e-06), make_tuple(-3, 2, 1, 1.1163526257063670e-06), make_tuple(-3, 2, 2, 1.1247885397545260e-07), make_tuple(-3, 2, 3, 9.3732378312877180e-10), make_tuple(-3, 3, -3, 7.8110315260730980e-12), make_tuple(-3, 3, -2, 9.3732378312877180e-10), make_tuple(-3, 3, -1, 9.3029385475530600e-09), make_tuple(-3, 3, 0, 1.8871452166992610e-08), make_tuple(-3, 3, 1, 9.3029385475530600e-09), make_tuple(-3, 3, 2, 9.3732378312877180e-10), make_tuple(-3, 3, 3, 7.8110315260730980e-12), make_tuple(-2, -3, -3, 9.3732378312877180e-10), make_tuple(-2, -3, -2, 1.1247885397545260e-07), make_tuple(-2, -3, -1, 1.1163526257063670e-06), make_tuple(-2, -3, 0, 2.2645742600391130e-06), make_tuple(-2, -3, 1, 1.1163526257063670e-06), make_tuple(-2, -3, 2, 1.1247885397545260e-07), make_tuple(-2, -3, 3, 9.3732378312877180e-10), make_tuple(-2, -2, -3, 1.1247885397545260e-07), make_tuple(-2, -2, -2, 1.3497462477054310e-05), make_tuple(-2, -2, -1, 1.3396231508476410e-04), make_tuple(-2, -2, 0, 2.7174891120469350e-04), make_tuple(-2, -2, 1, 1.3396231508476410e-04), make_tuple(-2, -2, 2, 1.3497462477054310e-05), make_tuple(-2, -2, 3, 1.1247885397545260e-07), make_tuple(-2, -1, -3, 1.1163526257063670e-06), make_tuple(-2, -1, -2, 1.3396231508476410e-04), make_tuple(-2, -1, -1, 1.3295759772162830e-03), make_tuple(-2, -1, 0, 2.6971079437065830e-03), make_tuple(-2, -1, 1, 1.3295759772162830e-03), make_tuple(-2, -1, 2, 1.3396231508476410e-04), make_tuple(-2, -1, 3, 1.1163526257063670e-06), make_tuple(-2, 0, -3, 2.2645742600391130e-06), make_tuple(-2, 0, -2, 2.7174891120469350e-04), make_tuple(-2, 0, -1, 2.6971079437065830e-03), make_tuple(-2, 0, 0, 5.4712114122544960e-03), make_tuple(-2, 0, 1, 2.6971079437065830e-03), make_tuple(-2, 0, 2, 2.7174891120469350e-04), make_tuple(-2, 0, 3, 2.2645742600391130e-06), make_tuple(-2, 1, -3, 1.1163526257063670e-06), make_tuple(-2, 1, -2, 1.3396231508476410e-04), make_tuple(-2, 1, -1, 1.3295759772162830e-03), make_tuple(-2, 1, 0, 2.6971079437065830e-03), make_tuple(-2, 1, 1, 1.3295759772162830e-03), make_tuple(-2, 1, 2, 1.3396231508476410e-04), make_tuple(-2, 1, 3, 1.1163526257063670e-06), make_tuple(-2, 2, -3, 1.1247885397545260e-07), make_tuple(-2, 2, -2, 1.3497462477054310e-05), make_tuple(-2, 2, -1, 1.3396231508476410e-04), make_tuple(-2, 2, 0, 2.7174891120469350e-04), make_tuple(-2, 2, 1, 1.3396231508476410e-04), make_tuple(-2, 2, 2, 1.3497462477054310e-05), make_tuple(-2, 2, 3, 1.1247885397545260e-07), make_tuple(-2, 3, -3, 9.3732378312877180e-10), make_tuple(-2, 3, -2, 1.1247885397545260e-07), make_tuple(-2, 3, -1, 1.1163526257063670e-06), make_tuple(-2, 3, 0, 2.2645742600391130e-06), make_tuple(-2, 3, 1, 1.1163526257063670e-06), make_tuple(-2, 3, 2, 1.1247885397545260e-07), make_tuple(-2, 3, 3, 9.3732378312877180e-10), make_tuple(-1, -3, -3, 9.3029385475530600e-09), make_tuple(-1, -3, -2, 1.1163526257063670e-06), make_tuple(-1, -3, -1, 1.1079799810135690e-05), make_tuple(-1, -3, 0, 2.2475899530888190e-05), make_tuple(-1, -3, 1, 1.1079799810135690e-05), make_tuple(-1, -3, 2, 1.1163526257063670e-06), make_tuple(-1, -3, 3, 9.3029385475530600e-09), make_tuple(-1, -2, -3, 1.1163526257063670e-06), make_tuple(-1, -2, -2, 1.3396231508476410e-04), make_tuple(-1, -2, -1, 1.3295759772162830e-03), make_tuple(-1, -2, 0, 2.6971079437065830e-03), make_tuple(-1, -2, 1, 1.3295759772162830e-03), make_tuple(-1, -2, 2, 1.3396231508476410e-04), make_tuple(-1, -2, 3, 1.1163526257063670e-06), make_tuple(-1, -1, -3, 1.1079799810135690e-05), make_tuple(-1, -1, -2, 1.3295759772162830e-03), make_tuple(-1, -1, -1, 1.3196041573871610e-02), make_tuple(-1, -1, 0, 2.6768796341287840e-02), make_tuple(-1, -1, 1, 1.3196041573871610e-02), make_tuple(-1, -1, 2, 1.3295759772162830e-03), make_tuple(-1, -1, 3, 1.1079799810135690e-05), make_tuple(-1, 0, -3, 2.2475899530888190e-05), make_tuple(-1, 0, -2, 2.6971079437065830e-03), make_tuple(-1, 0, -1, 2.6768796341287840e-02), make_tuple(-1, 0, 0, 5.4301773266625870e-02), make_tuple(-1, 0, 1, 2.6768796341287840e-02), make_tuple(-1, 0, 2, 2.6971079437065830e-03), make_tuple(-1, 0, 3, 2.2475899530888190e-05), make_tuple(-1, 1, -3, 1.1079799810135690e-05), make_tuple(-1, 1, -2, 1.3295759772162830e-03), make_tuple(-1, 1, -1, 1.3196041573871610e-02), make_tuple(-1, 1, 0, 2.6768796341287840e-02), make_tuple(-1, 1, 1, 1.3196041573871610e-02), make_tuple(-1, 1, 2, 1.3295759772162830e-03), make_tuple(-1, 1, 3, 1.1079799810135690e-05), make_tuple(-1, 2, -3, 1.1163526257063670e-06), make_tuple(-1, 2, -2, 1.3396231508476410e-04), make_tuple(-1, 2, -1, 1.3295759772162830e-03), make_tuple(-1, 2, 0, 2.6971079437065830e-03), make_tuple(-1, 2, 1, 1.3295759772162830e-03), make_tuple(-1, 2, 2, 1.3396231508476410e-04), make_tuple(-1, 2, 3, 1.1163526257063670e-06), make_tuple(-1, 3, -3, 9.3029385475530600e-09), make_tuple(-1, 3, -2, 1.1163526257063670e-06), make_tuple(-1, 3, -1, 1.1079799810135690e-05), make_tuple(-1, 3, 0, 2.2475899530888190e-05), make_tuple(-1, 3, 1, 1.1079799810135690e-05), make_tuple(-1, 3, 2, 1.1163526257063670e-06), make_tuple(-1, 3, 3, 9.3029385475530600e-09), make_tuple(0, -3, -3, 1.8871452166992610e-08), make_tuple(0, -3, -2, 2.2645742600391130e-06), make_tuple(0, -3, -1, 2.2475899530888190e-05), make_tuple(0, -3, 0, 4.5593428435454130e-05), make_tuple(0, -3, 1, 2.2475899530888190e-05), make_tuple(0, -3, 2, 2.2645742600391130e-06), make_tuple(0, -3, 3, 1.8871452166992610e-08), make_tuple(0, -2, -3, 2.2645742600391130e-06), make_tuple(0, -2, -2, 2.7174891120469350e-04), make_tuple(0, -2, -1, 2.6971079437065830e-03), make_tuple(0, -2, 0, 5.4712114122544960e-03), make_tuple(0, -2, 1, 2.6971079437065830e-03), make_tuple(0, -2, 2, 2.7174891120469350e-04), make_tuple(0, -2, 3, 2.2645742600391130e-06), make_tuple(0, -1, -3, 2.2475899530888190e-05), make_tuple(0, -1, -2, 2.6971079437065830e-03), make_tuple(0, -1, -1, 2.6768796341287840e-02), make_tuple(0, -1, 0, 5.4301773266625870e-02), make_tuple(0, -1, 1, 2.6768796341287840e-02), make_tuple(0, -1, 2, 2.6971079437065830e-03), make_tuple(0, -1, 3, 2.2475899530888190e-05), make_tuple(0, 0, -3, 4.5593428435454130e-05), make_tuple(0, 0, -2, 5.4712114122544960e-03), make_tuple(0, 0, -1, 5.4301773266625870e-02), make_tuple(0, 0, 0, 1.1015372310005720e-01), make_tuple(0, 0, 1, 5.4301773266625870e-02), make_tuple(0, 0, 2, 5.4712114122544960e-03), make_tuple(0, 0, 3, 4.5593428435454130e-05), make_tuple(0, 1, -3, 2.2475899530888190e-05), make_tuple(0, 1, -2, 2.6971079437065830e-03), make_tuple(0, 1, -1, 2.6768796341287840e-02), make_tuple(0, 1, 0, 5.4301773266625870e-02), make_tuple(0, 1, 1, 2.6768796341287840e-02), make_tuple(0, 1, 2, 2.6971079437065830e-03), make_tuple(0, 1, 3, 2.2475899530888190e-05), make_tuple(0, 2, -3, 2.2645742600391130e-06), make_tuple(0, 2, -2, 2.7174891120469350e-04), make_tuple(0, 2, -1, 2.6971079437065830e-03), make_tuple(0, 2, 0, 5.4712114122544960e-03), make_tuple(0, 2, 1, 2.6971079437065830e-03), make_tuple(0, 2, 2, 2.7174891120469350e-04), make_tuple(0, 2, 3, 2.2645742600391130e-06), make_tuple(0, 3, -3, 1.8871452166992610e-08), make_tuple(0, 3, -2, 2.2645742600391130e-06), make_tuple(0, 3, -1, 2.2475899530888190e-05), make_tuple(0, 3, 0, 4.5593428435454130e-05), make_tuple(0, 3, 1, 2.2475899530888190e-05), make_tuple(0, 3, 2, 2.2645742600391130e-06), make_tuple(0, 3, 3, 1.8871452166992610e-08), make_tuple(1, -3, -3, 9.3029385475530600e-09), make_tuple(1, -3, -2, 1.1163526257063670e-06), make_tuple(1, -3, -1, 1.1079799810135690e-05), make_tuple(1, -3, 0, 2.2475899530888190e-05), make_tuple(1, -3, 1, 1.1079799810135690e-05), make_tuple(1, -3, 2, 1.1163526257063670e-06), make_tuple(1, -3, 3, 9.3029385475530600e-09), make_tuple(1, -2, -3, 1.1163526257063670e-06), make_tuple(1, -2, -2, 1.3396231508476410e-04), make_tuple(1, -2, -1, 1.3295759772162830e-03), make_tuple(1, -2, 0, 2.6971079437065830e-03), make_tuple(1, -2, 1, 1.3295759772162830e-03), make_tuple(1, -2, 2, 1.3396231508476410e-04), make_tuple(1, -2, 3, 1.1163526257063670e-06), make_tuple(1, -1, -3, 1.1079799810135690e-05), make_tuple(1, -1, -2, 1.3295759772162830e-03), make_tuple(1, -1, -1, 1.3196041573871610e-02), make_tuple(1, -1, 0, 2.6768796341287840e-02), make_tuple(1, -1, 1, 1.3196041573871610e-02), make_tuple(1, -1, 2, 1.3295759772162830e-03), make_tuple(1, -1, 3, 1.1079799810135690e-05), make_tuple(1, 0, -3, 2.2475899530888190e-05), make_tuple(1, 0, -2, 2.6971079437065830e-03), make_tuple(1, 0, -1, 2.6768796341287840e-02), make_tuple(1, 0, 0, 5.4301773266625870e-02), make_tuple(1, 0, 1, 2.6768796341287840e-02), make_tuple(1, 0, 2, 2.6971079437065830e-03), make_tuple(1, 0, 3, 2.2475899530888190e-05), make_tuple(1, 1, -3, 1.1079799810135690e-05), make_tuple(1, 1, -2, 1.3295759772162830e-03), make_tuple(1, 1, -1, 1.3196041573871610e-02), make_tuple(1, 1, 0, 2.6768796341287840e-02), make_tuple(1, 1, 1, 1.3196041573871610e-02), make_tuple(1, 1, 2, 1.3295759772162830e-03), make_tuple(1, 1, 3, 1.1079799810135690e-05), make_tuple(1, 2, -3, 1.1163526257063670e-06), make_tuple(1, 2, -2, 1.3396231508476410e-04), make_tuple(1, 2, -1, 1.3295759772162830e-03), make_tuple(1, 2, 0, 2.6971079437065830e-03), make_tuple(1, 2, 1, 1.3295759772162830e-03), make_tuple(1, 2, 2, 1.3396231508476410e-04), make_tuple(1, 2, 3, 1.1163526257063670e-06), make_tuple(1, 3, -3, 9.3029385475530600e-09), make_tuple(1, 3, -2, 1.1163526257063670e-06), make_tuple(1, 3, -1, 1.1079799810135690e-05), make_tuple(1, 3, 0, 2.2475899530888190e-05), make_tuple(1, 3, 1, 1.1079799810135690e-05), make_tuple(1, 3, 2, 1.1163526257063670e-06), make_tuple(1, 3, 3, 9.3029385475530600e-09), make_tuple(2, -3, -3, 9.3732378312877180e-10), make_tuple(2, -3, -2, 1.1247885397545260e-07), make_tuple(2, -3, -1, 1.1163526257063670e-06), make_tuple(2, -3, 0, 2.2645742600391130e-06), make_tuple(2, -3, 1, 1.1163526257063670e-06), make_tuple(2, -3, 2, 1.1247885397545260e-07), make_tuple(2, -3, 3, 9.3732378312877180e-10), make_tuple(2, -2, -3, 1.1247885397545260e-07), make_tuple(2, -2, -2, 1.3497462477054310e-05), make_tuple(2, -2, -1, 1.3396231508476410e-04), make_tuple(2, -2, 0, 2.7174891120469350e-04), make_tuple(2, -2, 1, 1.3396231508476410e-04), make_tuple(2, -2, 2, 1.3497462477054310e-05), make_tuple(2, -2, 3, 1.1247885397545260e-07), make_tuple(2, -1, -3, 1.1163526257063670e-06), make_tuple(2, -1, -2, 1.3396231508476410e-04), make_tuple(2, -1, -1, 1.3295759772162830e-03), make_tuple(2, -1, 0, 2.6971079437065830e-03), make_tuple(2, -1, 1, 1.3295759772162830e-03), make_tuple(2, -1, 2, 1.3396231508476410e-04), make_tuple(2, -1, 3, 1.1163526257063670e-06), make_tuple(2, 0, -3, 2.2645742600391130e-06), make_tuple(2, 0, -2, 2.7174891120469350e-04), make_tuple(2, 0, -1, 2.6971079437065830e-03), make_tuple(2, 0, 0, 5.4712114122544960e-03), make_tuple(2, 0, 1, 2.6971079437065830e-03), make_tuple(2, 0, 2, 2.7174891120469350e-04), make_tuple(2, 0, 3, 2.2645742600391130e-06), make_tuple(2, 1, -3, 1.1163526257063670e-06), make_tuple(2, 1, -2, 1.3396231508476410e-04), make_tuple(2, 1, -1, 1.3295759772162830e-03), make_tuple(2, 1, 0, 2.6971079437065830e-03), make_tuple(2, 1, 1, 1.3295759772162830e-03), make_tuple(2, 1, 2, 1.3396231508476410e-04), make_tuple(2, 1, 3, 1.1163526257063670e-06), make_tuple(2, 2, -3, 1.1247885397545260e-07), make_tuple(2, 2, -2, 1.3497462477054310e-05), make_tuple(2, 2, -1, 1.3396231508476410e-04), make_tuple(2, 2, 0, 2.7174891120469350e-04), make_tuple(2, 2, 1, 1.3396231508476410e-04), make_tuple(2, 2, 2, 1.3497462477054310e-05), make_tuple(2, 2, 3, 1.1247885397545260e-07), make_tuple(2, 3, -3, 9.3732378312877180e-10), make_tuple(2, 3, -2, 1.1247885397545260e-07), make_tuple(2, 3, -1, 1.1163526257063670e-06), make_tuple(2, 3, 0, 2.2645742600391130e-06), make_tuple(2, 3, 1, 1.1163526257063670e-06), make_tuple(2, 3, 2, 1.1247885397545260e-07), make_tuple(2, 3, 3, 9.3732378312877180e-10), make_tuple(3, -3, -3, 7.8110315260730980e-12), make_tuple(3, -3, -2, 9.3732378312877180e-10), make_tuple(3, -3, -1, 9.3029385475530600e-09), make_tuple(3, -3, 0, 1.8871452166992610e-08), make_tuple(3, -3, 1, 9.3029385475530600e-09), make_tuple(3, -3, 2, 9.3732378312877180e-10), make_tuple(3, -3, 3, 7.8110315260730980e-12), make_tuple(3, -2, -3, 9.3732378312877180e-10), make_tuple(3, -2, -2, 1.1247885397545260e-07), make_tuple(3, -2, -1, 1.1163526257063670e-06), make_tuple(3, -2, 0, 2.2645742600391130e-06), make_tuple(3, -2, 1, 1.1163526257063670e-06), make_tuple(3, -2, 2, 1.1247885397545260e-07), make_tuple(3, -2, 3, 9.3732378312877180e-10), make_tuple(3, -1, -3, 9.3029385475530600e-09), make_tuple(3, -1, -2, 1.1163526257063670e-06), make_tuple(3, -1, -1, 1.1079799810135690e-05), make_tuple(3, -1, 0, 2.2475899530888190e-05), make_tuple(3, -1, 1, 1.1079799810135690e-05), make_tuple(3, -1, 2, 1.1163526257063670e-06), make_tuple(3, -1, 3, 9.3029385475530600e-09), make_tuple(3, 0, -3, 1.8871452166992610e-08), make_tuple(3, 0, -2, 2.2645742600391130e-06), make_tuple(3, 0, -1, 2.2475899530888190e-05), make_tuple(3, 0, 0, 4.5593428435454130e-05), make_tuple(3, 0, 1, 2.2475899530888190e-05), make_tuple(3, 0, 2, 2.2645742600391130e-06), make_tuple(3, 0, 3, 1.8871452166992610e-08), make_tuple(3, 1, -3, 9.3029385475530600e-09), make_tuple(3, 1, -2, 1.1163526257063670e-06), make_tuple(3, 1, -1, 1.1079799810135690e-05), make_tuple(3, 1, 0, 2.2475899530888190e-05), make_tuple(3, 1, 1, 1.1079799810135690e-05), make_tuple(3, 1, 2, 1.1163526257063670e-06), make_tuple(3, 1, 3, 9.3029385475530600e-09), make_tuple(3, 2, -3, 9.3732378312877180e-10), make_tuple(3, 2, -2, 1.1247885397545260e-07), make_tuple(3, 2, -1, 1.1163526257063670e-06), make_tuple(3, 2, 0, 2.2645742600391130e-06), make_tuple(3, 2, 1, 1.1163526257063670e-06), make_tuple(3, 2, 2, 1.1247885397545260e-07), make_tuple(3, 2, 3, 9.3732378312877180e-10), make_tuple(3, 3, -3, 7.8110315260730980e-12), make_tuple(3, 3, -2, 9.3732378312877180e-10), make_tuple(3, 3, -1, 9.3029385475530600e-09), make_tuple(3, 3, 0, 1.8871452166992610e-08), make_tuple(3, 3, 1, 9.3029385475530600e-09), make_tuple(3, 3, 2, 9.3732378312877180e-10), make_tuple(3, 3, 3, 7.8110315260730980e-12)};
	}

	inline static const vector3<O> approximateGradient(const int &i, const int &j, const int &k, const shift_invariant_space3< tp3cubic, O, I> *L) {
		O sx = O(1./(L->getScale()*12.));

		O dx = (
		 -1.*L->GV(i-2, j, k) + 
		  8.*L->GV(i-1, j, k) +
		 -8.*L->GV(i+1, j, k) + 
			 L->GV(i+2, j, k)) * sx;
		
		O dy = (
		 (-1.*L->GV(i, j-2, k)) + 
		 ( 8.*L->GV(i, j-1, k)) + 
		 (-8.*L->GV(i, j+1, k)) + 
		 (	L->GV(i, j+2, k))) * sx;

		O dz = (
		 -1.*L->GV(i, j, k-2) + 
		  8.*L->GV(i, j, k-1) + 
		 -8.*L->GV(i, j, k+1) + 
			 L->GV(i, j, k+2)) * sx; 
		return vector3<O>(dx,dy,dz);
	}

	// This does the actual semi-descrete convolution sum
	// the idea behind this, is to allow basis functions
	// to provide a potentially optimized version of
	// the convolution sum for particular lattices.
	inline static const O convolutionSum(const vector3<I> &p, const shift_invariant_space3< tp3cubic, O, I> *lattice) {
		O sum = 0;

		vector3<I> vox;
		
		I dh = lattice->getScale();

		vox.i = p.i/dh;
		vox.j = p.j/dh;
		vox.k = p.k/dh;

		int 
			vx = (int)floor(vox.i), 
			vy = (int)floor(vox.j), 
			vz = (int)floor(vox.k);

		for(int i = -1; i <= 2; i++)
			for(int j = -1; j <= 2; j++)
				for(int k = -1; k <= 2; k++){
					sum += lattice->GV(vx + i, vy + j, vz + k) *
							tp3cubic::M(vox.i - I(vx + i), 
								vox.j - I(vy + j),
								vox.k - I(vz + k));
				}
		return sum;		
	}
	inline static const vector3<O> convolutionSumNormal(const vector3<I> &p, shift_invariant_space3< tp3cubic, O, I> *lattice) {
		vector3<O> sum = vector3<O>(0,0,0);

		vector3<I> vox;
		
		I dh = lattice->getScale();

		vox.i = p.i/dh;
		vox.j = p.j/dh;
		vox.k = p.k/dh;

		int 
			vx = (int)floor(vox.i), 
			vy = (int)floor(vox.j), 
			vz = (int)floor(vox.k);

		for(int i = -1; i <= 2; i++)
			for(int j = -1; j <= 2; j++)
				for(int k = -1; k <= 2; k++){
					sum = sum + lattice->GN(vx + i, vy + j, vz + k) *
							tp3cubic::M(vox.i - I(vx + i), 
								vox.j - I(vy + j),
								vox.k - I(vz + k));
				}
		return sum;		
	}
	static O interpolationFilter(const I &u, const I &v, const I &w, const I &h) {
		return 1./((1./27.)*(2.+cos(u))*(2.+cos(v))*(2.+cos(w)));
	}

	static O poissonFilter(const I &u, const I &v, const I &w, const I &h) {
		O invh = 1./(h*h);
		O lambda = (2./3.)*
							( invh*(-7. + cos(u))*sin(0.5*u)*sin(0.5*u) +
								invh*(-7. + cos(v))*sin(0.5*v)*sin(0.5*v) +
								invh*(-7. + cos(w))*sin(0.5*w)*sin(0.5*w) );
		return 	1./lambda;
	}
	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		divergenceFilter(L *x, L *y, L *z){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			double sx = -1./(x->getScale()*12.);

			double dx = (
			 -1.*x->GV(i-2, j, k) + 
			  8.*x->GV(i-1, j, k) +
			 -8.*x->GV(i+1, j, k) + 
				 x->GV(i+2, j, k)) * sx;
		
			double dy = (
			 (-1.*y->GV(i, j-2, k)) + 
			 ( 8.*y->GV(i, j-1, k)) + 
			 (-8.*y->GV(i, j+1, k)) + 
			 (	y->GV(i, j+2, k))) * sx;

			double dz = (
			 -1.*z->GV(i, j, k-2) + 
			  8.*z->GV(i, j, k-1) + 
			 -8.*z->GV(i, j, k+1) + 
				 z->GV(i, j, k+2)) * sx; 

			return dx+dy+dz;
		};
	}

	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		shiftedDivergenceFilter(L *x, L *y, L *z){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			double sx = -1./(x->getScale()*24.);

			double dx = (
				 -1.*x->GV(i-1, j, k) + 
				  27.*x->GV(i, j, k) +
				 -27.*x->GV(i+1, j, k) + 
					 x->GV(i+2, j, k)) * sx;
					
			double dy = (
				 (-1.*y->GV(i, j-1, k)) + 
				 ( 27.*y->GV(i, j, k)) + 
				 (-27.*y->GV(i, j+1, k)) + 
				 (	y->GV(i, j+2, k))) * sx;

			double dz = (
				 -1.*z->GV(i, j, k-1) + 
				  27.*z->GV(i, j, k) + 
				 -27.*z->GV(i, j, k+1) + 
					 z->GV(i, j, k+2)) * sx; 

			return dx+dy+dz;
		};
	}

	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		dxFilter0bd(L *x){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			O sx = -1./(x->getScale()*12.);
			return (
			 -1.*x->GV(i-2, j, k) + 
			  8.*x->GV(i-1, j, k) +
			 -8.*x->GV(i+1, j, k) + 
				 x->GV(i+2, j, k)) * sx;

		};
	}

	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		dyFilter0bd(L *y){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			O sy = -1./(y->getScale()*12.);
			return (
			 (-1.*y->GV(i, j-2, k)) + 
			 ( 8.*y->GV(i, j-1, k)) + 
			 (-8.*y->GV(i, j+1, k)) + 
			 (	y->GV(i, j+2, k))) * sy;
		};
	}

	template <class L>
	static std::function<O(const int &, const int &, const int &, std::function<O(const int &, const int &, const int & )>)> 
		dzFilter0bd(L *z){
		return [&](const int & i, const int &j, const int &k, std::function<O(const int &, const int &, const int & )> s){
			O sz = -1./(z->getScale()*12.);
			return (
			 (-1.*z->GV(i, j-2, k)) + 
			 ( 8.*z->GV(i, j-1, k)) + 
			 (-8.*z->GV(i, j+1, k)) + 
			 (	z->GV(i, j+2, k))) * sz;
		};
	}
};
};
#endif //_TP_BSPLINE3_H_