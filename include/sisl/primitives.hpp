/**
 * Some additional ``primitive'' data types.
 * 
 * @author Joshua Horacsek
 */

#include <cmath>

#ifndef _PRIMITIVES_H
#define _PRIMITIVES_H

namespace sisl{
template <class T = float>
struct vector3{
	T i,j,k;

	inline vector3(){i=j=k=0.;}
	inline vector3(T ii, T jj, T kk){i=ii;j=jj;k=kk;}

	inline vector3<T> &operator+=(const vector3<T> &v) {
		this->i += v.i; this->j += v.j; this->k += v.k; 
		return *this;
	}
	inline vector3<T> &operator-=(const vector3<T> &v) {
		this->i -= v.i; this->j -= v.j; this->k -= v.k; 
		return *this;
	}
	inline vector3<T> operator +(const vector3<T> &v) const {
		return vector3<T>(i + v.i, j + v.j, k + v.k);
	}
	inline vector3<T> operator -(const vector3<T> &v) const {
		return vector3<T>(i - v.i, j - v.j, k - v.k);
	}
	inline T operator*(const vector3 &v) const {
		return i*v.i + j*v.j + k * v.k;
	}
	inline vector3<T> operator%(const vector3<T> &v) const {
		return vector3<T>(j*v.k - k*v.j, v.i*k - i*v.k, i*v.j - j*v.i);
	}
	inline vector3<T> operator -() const{
		return vector3<T>(-i, -j, -k);
	}
	inline vector3<T> operator*(T scalar) const{
		return vector3 (i * scalar, j * scalar, k * scalar);
	}
	inline vector3<T> & normalize() {
		T norm2 = ((T)1.0)/((T)sqrt((*this) * (*this)));
		i *= norm2;
		j *= norm2;
		k *= norm2;
		return *this;
	}

	inline T length2() const{
		return *this * *this;
	}

	inline T norm() const {
		return ((T)sqrt((*this) * (*this)));
	}
};

// Column, Row indexed, not sure why to be honest...
template <class T = float>
struct matrix4x4 {
	T mat[4][4];
	matrix4x4(){this->identity();}

	void identity(){
		mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;
		mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = 0;
		mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1; mat[2][3] = 0;
		mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;
	}

	void zero() {
		mat[0][0] = 0; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;
		mat[1][0] = 0; mat[1][1] = 0; mat[1][2] = 0; mat[1][3] = 0;
		mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 0; mat[2][3] = 0;
		mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 0;
	}

	vector3<T> operator*(const vector3<T> &v) const {
		vector3<T> r;
		T res; 
		r.i = mat[0][0] * v.i + mat[1][0] * v.j + mat[2][0] * v.k + mat[3][0];
		r.j = mat[0][1] * v.i + mat[1][1] * v.j + mat[2][1] * v.k + mat[3][1];
		r.k = mat[0][2] * v.i + mat[1][2] * v.j + mat[2][2] * v.k + mat[3][2];
		res = mat[0][3] * v.i + mat[1][3] * v.j + mat[2][3] * v.k + mat[3][3];
		return r* (1./res);
	}

	matrix4x4<T> operator*(const matrix4x4<T> &m) const{
		matrix4x4<T> n;
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++) {
				n.mat[j][i] = 0;
				for(int k = 0; k < 4; k++)
					n.mat[j][i] += mat[k][i]*m.mat[j][k];
			}

		return n;
	}

	matrix4x4<T> operator*(const T scalar) const{
		matrix4x4<T> n;
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++) {
				n.mat[j][i] = mat[j][i] * scalar;
			}

		return n;
	}

	matrix4x4<T> operator+(const matrix4x4<T> &m) const{
		matrix4x4<T> n;
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++) {
				n.mat[j][i] = mat[j][i] + m.mat[j][i];
			}

		return n;
	}

	matrix4x4<T> &operator+=(const matrix4x4<T> &m) {
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++) {
				mat[j][i] += m.mat[j][i];
			}
		return *this;
	}

	matrix4x4<T> &operator-=(const matrix4x4<T> &m) {
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++) {
				mat[j][i] -= m.mat[j][i];
			}
		return *this;
	}

	static matrix4x4<T> translate(T x, T y, T z){
		matrix4x4<T> r;
		r.mat[3][0] = x;
		r.mat[3][1] = y;
		r.mat[3][2] = z;
		return r;
	}

	static matrix4x4<T> translate(vector3<T> p){
		matrix4x4<T> r;
		r.mat[3][0] = p.i;
		r.mat[3][1] = p.j;
		r.mat[3][2] = p.k;
		return r;
	}

	static matrix4x4<T> scale(T x, T y, T z){
		matrix4x4<T> r;
		r.mat[0][0] = x;
		r.mat[1][1] = y;
		r.mat[2][2] = z;
		return r;
	}

	static matrix4x4<T> uniformScale(T x){
		matrix4x4<T> r;
		r.mat[0][0] = x;
		r.mat[1][1] = x;
		r.mat[2][2] = x;
		return r;
	}
};

template <class T = float>
struct vertex3 {
	vector3<T> p;
	vector3<T> n;

	vertex3(vector3<T> _p, vector3<T> _n) : p(_p), n(_n) {}
	vertex3(T x, T y, T z, T i, T j, T k) :  p(vector3<T>(x,y,z)), n(vector3<T>(i,j,k)) {} 
};

template <class T = float>
struct triangle {
	vertex3<T> v[3];
};

};

#endif // _PRIMITIVES_H
