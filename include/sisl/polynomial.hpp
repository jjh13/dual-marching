#ifndef __NVAR_POLYNOMIAL__
#define __NVAR_POLYNOMIAL__

#include <vector>
#include <cmath>
#include <stdexcept>

namespace sisl {

	// Forward Declarations
	template<int N, class I, class O> class Polynomial;

	// Base class for a term in a polynomial.
	// The design is simple, Constants are terms,
	// a variable is a term, and a polynomial is a term
	// This lends itself easily to a recursive evaluation
	// scheme.
	template<int N, class I=double, class O=double>
	class Term {
	public:
		virtual ~Term() { }
		virtual O evaluate(const std::vector<I> &) const = 0;
		virtual void expand() = 0;
		virtual void collect_like_terms() = 0;
		virtual Term<N,I,O> *diff(const int &) = 0;
		virtual bool is_const() = 0;
		virtual bool is_const(const int &) = 0;
		virtual Term<N,I,O> *clone() const = 0;
	};


	template<int N, class I=double, class O=double>
	class ConstantTerm: public Term<N,I,O> {
	public:
		ConstantTerm(const O &c) : constant(c) {}
		ConstantTerm(const ConstantTerm<N,I,O> &c){ this->constant = c.constant; }
		virtual ~ConstantTerm(){ }

		ConstantTerm<N,I,O> &operator=(const ConstantTerm <N,I,O> &rhs) {
			this->constant = rhs.constant;
			return *this;
		}

		virtual O evaluate(const std::vector<I> &) const { return this->constant; }
		virtual void expand() { }
		virtual void collect_like_terms() { };
		virtual Term<N,I,O> *diff(const int &) { return new ConstantTerm<N,I,O>((O)0); }
		virtual bool is_const() { return true; }
		virtual bool is_const(const int &) { return true; }
		virtual ConstantTerm<N,I,O> *clone() const { return new ConstantTerm<N,I,O>(*this); }
	private:
		O constant;
	};

	template<int N, class I=double, class O=double>
	class VariableTerm: public Term<N,I,O> {
	public: 
		VariableTerm(const int &termIndex) {
			if(termIndex >= N) throw new std::range_error("Variable exceeded polynomial dimensionality!");
			this->variable = termIndex;
		}
		VariableTerm(const VariableTerm<N,I,O> &vt) {
			if(vt.variable >= N) throw new std::range_error("Variable exceeded polynomial dimensionality!");
			this->variable = vt.variable;
		}
		virtual ~VariableTerm(){ }

		VariableTerm<N,I,O> &operator=(const VariableTerm <N,I,O> &rhs) {
			this->variable = rhs.variable;
			return *this;
		}

		virtual O evaluate(const std::vector<I> &x) const { return x[this->variable]; }
		virtual void expand() {	}
		virtual void collect_like_terms() { };
		virtual Term<N,I,O> *diff(const int &di) { return new ConstantTerm<N,I,O>(di==this->variable?(O)1:(O)0); }
		virtual bool is_const() { return false; }
		virtual bool is_const(const int &di) { return di==this->variable; }
		virtual VariableTerm<N,I,O> *clone() const {
			return new VariableTerm<N,I,O>(*this);
		}
	private:
		int variable;
	};

	// Small ancilary class for declaring variables 
	// in expressions.
	// I.E. Polynomial<3> p = 100*PVar<3>(0)*PVar<3>(1) 
	template<int N, class I=double, class O=double>
	class PVar {
	public:
		PVar(int v) : variable(v){ }
		int getVariable() const {return variable;}
		operator const Polynomial<N,I,O>() const {
			return Polynomial<N,I,O>(*this);
		} 
	private:
		int variable;
	};

	// Class PVVec

	template<int N, class I=double, class O=double>
	class Polynomial: public Term<N,I,O> {
	public:
		Polynomial() : Polynomial((O)0) { }
		Polynomial(O i) : dirty(false), expanded(true) {
			poly_term t = { new ConstantTerm<N,I,O>(i), 1 };
			std::vector<poly_term> tlist;
			this->terms = new std::vector<std::vector<poly_term> >();
			tlist.push_back(t);
			this->terms->push_back(tlist);
		}

		Polynomial(const PVar<N,I,O> &p) : dirty(false), expanded(true) {
			poly_term t = { new VariableTerm<N,I,O>(p.getVariable()), 1 };
			std::vector<poly_term> tlist;
			this->terms = new std::vector<std::vector<poly_term> >();
			tlist.push_back(t);
			this->terms->push_back(tlist);
		}

		Polynomial(const Polynomial<N,I,O> &p) : dirty(false), expanded(true) {
			this->terms = new std::vector<std::vector<poly_term> >();
			for(int i = 0; i < (*p.terms).size(); i++) {
				std::vector<poly_term> tlist;
				for(int j = 0; j < (*p.terms)[i].size(); j++) {
					poly_term t_prime = {(*p.terms)[i][j].term->clone(), (*p.terms)[i][j].power};
					tlist.push_back(t_prime);
				}
				this->terms->push_back(tlist);
			}
			this->dirty = p.dirty;
			this->expanded = p.expanded;
		}

		Polynomial(Term<N,I,O> *trm) : dirty(false), expanded(true) {
			poly_term t = { trm, 1 };
			std::vector<poly_term> tlist;
			this->terms = new std::vector<std::vector<poly_term> >();
			tlist.push_back(t);
			this->terms->push_back(tlist);
		}

		~Polynomial() {
			for(int i = 0; i < (*terms).size(); i++)
				for(int j = 0; j < (*terms)[i].size(); j++)
					delete (*terms)[i][j].term;
			delete this->terms;
		}

		virtual O evaluate(const std::vector<I> &x) const {
			O accum = (O)0;
			for(int i = 0; i < terms->size(); i++) {
				O paccum = (O)1;
				for(int j = 0; j < (*terms)[i].size(); j++)
				 	paccum *= (O)std::pow((O)(*terms)[i][j].term->evaluate(x), (int)(*terms)[i][j].power);
				accum += paccum;
			}
			return accum;
		}

		virtual Term<N,I,O> *diff(const int &) { return this; }
		virtual bool is_const() { return false; };
		virtual bool is_const(const int &) { return false; }
		virtual Polynomial<N,I,O> *clone() const {
			return new Polynomial<N,I,O>(*this);
		}

		Polynomial<N,I,O> &operator=(Polynomial <N,I,O> rhs) {
			std::vector<std::vector<poly_term> > *tterms = this->terms;
			this->terms = rhs.terms;
			rhs.terms = tterms;
			this->dirty = rhs.dirty;
			this->expanded = rhs.expanded;
			return *this;
		}

		Polynomial<N,I,O> &operator+=(const Polynomial<N,I,O> &rhs) {
			for(int i = 0; i < (*rhs.terms).size(); i++) {
				std::vector<poly_term> tlist;
				for(int j = 0; j < (*rhs.terms)[i].size(); j++) {
					poly_term t_prime = {(*rhs.terms)[i][j].term->clone(), (*rhs.terms)[i][j].power};
					tlist.push_back(t_prime);
				}
				this->terms->push_back(tlist);
			}

			this->dirty = true;
			this->expanded = false;
			return *this;
		}

		Polynomial<N,I,O> &operator*=(const Polynomial<N,I,O> &rhs) {
			Polynomial<N,I,O> r = (O)1;

			poly_term p = {this->clone(), 1};
			poly_term q = {rhs.clone(), 1};

			(*r.terms)[0].push_back(p);
			(*r.terms)[0].push_back(q);

			*this = r;
			
			return *this;
		}

		friend inline Polynomial<N,I,O> operator+(Polynomial<N,I,O> lhs, const Polynomial<N,I,O> &rhs) {
			lhs += rhs;
			return lhs;
		}

		friend inline Polynomial<N,I,O> operator*(Polynomial<N,I,O> lhs, const Polynomial<N,I,O> &rhs) {
			lhs *= rhs;
			return lhs;
		}

		friend inline Polynomial<N,I,O> operator+(Polynomial<N,I,O> lhs, const PVar<N,I,O> &rhs) { return lhs + Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator+(const PVar<N,I,O> &rhs, Polynomial<N,I,O> lhs) { return lhs + Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator+(Polynomial<N,I,O> lhs, const O &rhs) { return lhs + Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator+(const O &rhs, Polynomial<N,I,O> lhs) { return lhs + Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator+(O lhs, const PVar<N,I,O> &rhs) { return Polynomial<N,I,O>(lhs) + Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator+(PVar<N,I,O> lhs, const O &rhs) { return Polynomial<N,I,O>(lhs) + Polynomial<N,I,O>(rhs); }

		friend inline Polynomial<N,I,O> operator*(Polynomial<N,I,O> lhs, const PVar<N,I,O> &rhs) { return lhs * Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator*(const PVar<N,I,O> &rhs, Polynomial<N,I,O> lhs) { return lhs * Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator*(Polynomial<N,I,O> lhs, const O &rhs) { return lhs * Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator*(const O &rhs, Polynomial<N,I,O> lhs) { return lhs * Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator*(O lhs, const PVar<N,I,O> &rhs) { return Polynomial<N,I,O>(lhs) * Polynomial<N,I,O>(rhs); }
		friend inline Polynomial<N,I,O> operator*(PVar<N,I,O> lhs, const O &rhs) { return Polynomial<N,I,O>(lhs) * Polynomial<N,I,O>(rhs); }


		virtual void expand() {

		}
		
		virtual void collect_like_terms() {

		}

	private:
		struct poly_term {
			Term<N,I,O> *term;
			unsigned int power;
		};

		std::vector<std::vector<poly_term> > *terms;
		bool dirty;
		bool expanded;
	};

	template<int N, class I, class O> Polynomial<N,I,O> operator +(const O &lhs, const PVar<N,I,O> &rhs) { return Polynomial<N,I,O>(lhs) + Polynomial<N,I,O>(rhs); }
	template<int N, class I, class O> Polynomial<N,I,O> operator +(const PVar<N,I,O> &lhs, const O &rhs) { return Polynomial<N,I,O>(lhs) + Polynomial<N,I,O>(rhs); }
	template<int N, class I, class O> Polynomial<N,I,O> operator +(const O &lhs, const O &rhs) { return Polynomial<N,I,O>(lhs) + Polynomial<N,I,O>(rhs); }

	template<int N, class I, class O> Polynomial<N,I,O> operator *(const O &lhs, const PVar<N,I,O> &rhs) { return Polynomial<N,I,O>(lhs) * Polynomial<N,I,O>(rhs); }
	template<int N, class I, class O> Polynomial<N,I,O> operator *(const PVar<N,I,O> &lhs, const O &rhs) { return Polynomial<N,I,O>(lhs) * Polynomial<N,I,O>(rhs); }
	template<int N, class I, class O> Polynomial<N,I,O> operator *(const PVar<N,I,O> &lhs, const PVar<N,I,O> &rhs) { return Polynomial<N,I,O>(lhs) * Polynomial<N,I,O>(rhs); }
}

#endif // __NVAR_POLYNOMIAL__