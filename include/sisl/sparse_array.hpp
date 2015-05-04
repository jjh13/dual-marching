/**
 * Generic class for handling 3-dimensional arrays, automatically
 * aligned for FFTW
 *
 * @author Joshua Horacsek
 *
 **/
#ifndef _SPARSE_ARRAY3_H
#define _SPARSE_ARRAY3_H

#include <string>
#include <memory>
#include <unordered_map>
#include <tuple>

#include <sisl/array.hpp>

// Useful for debugging
//#define NO_ARRAY_CHECK
namespace sisl {

template <class T, class Allocator = std::allocator<T> >
class sparse_array4 : array4<T,Allocator> {

private:
	T *_array;
	unsigned int _nx, _ny, _nz, _nw;
	std::unordered_map<int, T> siteMap;
	T defaultValue;

public:
	sparse_array4(const T &dv) : _nx(0), _ny(0), _nz(0), _nw(0), defaultValue(dv) {
		siteMap = std::unordered_map<int, T>();
	}
	sparse_array4(const int &nx, const int &ny, const int &nz, const int &nw, const T &dv) :  _nx(nx), _ny(ny), _nz(nz), _nw(nw), defaultValue(dv) { }
	~sparse_array4() {	}

	inline sparse_array4& operator=(const sparse_array4 &rhs) {
		throw "Not implemented";
	}

	inline void getDims(unsigned int *nx, unsigned int *ny, unsigned int *nz, unsigned int *nw) const {
		if(nx != NULL)
			*nx = _nx;
		if(ny != NULL)
			*ny = _ny;
		if(nz != NULL)
			*nz = _nz;
		if(nw != NULL)
			*nw = _nw;
	}

	inline T& operator()(int x, int y, int z, int w) {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
		if(!(z >= 0 && z < _nz)) throw "Index Z out of bounds!";
		if(!(w >= 0 && w < _nw)) throw "Index W out of bounds!";
#endif 
		int index = lIndex(x,y,z,w);
		auto iter = siteMap.find(index);
		if(iter == siteMap.end())
			siteMap[index] = defaultValue;

		return siteMap[index];
	} 

	inline const T& operator()(int x, int y, int z, int w) const {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
		if(!(z >= 0 && z < _nz)) throw "Index Z out of bounds!";
		if(!(w >= 0 && w < _nw)) throw "Index W out of bounds!";
#endif 
		int index = lIndex(x,y,z,w);
		auto iter = siteMap.find(index);
		if(iter == siteMap.end())
			return defaultValue;
		return iter->second;
	}

	inline bool hasValue(const int &x, const int &y, const int &z, const int &w){
		return siteMap.find(lIndex(x,y,z,w)) != siteMap.end();
	}

	inline int lIndex(const int &x, const int &y, const int &z, const int &w) const {
		return (z + _nz * (y + _ny*x))*_nw + w;
	}

	/**
	 * Gets a pointer to the memory address
	 * that the actual array is located at.
	 * FFTW will need this eventually.
	 **/
	T* getArray() const {
		return _array;
	}

	int size() const {
		return _nx*_ny*_nz*_nw;
	}

};

template <class T, class Allocator = std::allocator<T> >
class sparse_array3 : array3<T,Allocator> {

private:
	T *_array;
	unsigned int _nx, _ny, _nz;
	std::unordered_map<int, T> siteMap;
	T defaultValue;

public:
	sparse_array3(const T &dv) : _nx(0), _ny(0), _nz(0), defaultValue(dv) {
		siteMap = std::unordered_map<int, T>(); 
	}
	sparse_array3(const int &nx, const int &ny, const int &nz, const T &dv) :  _nx(nx), _ny(ny), _nz(nz), defaultValue(dv) { }
	~sparse_array3() { }

	inline sparse_array3& operator=(const sparse_array3 &rhs) {
		throw "Not implemented";
	}

	inline void getDims(unsigned int *nx, unsigned int *ny, unsigned int *nz) const {
		if(nx != NULL)
			*nx = _nx;
		if(ny != NULL)
			*ny = _ny;
		if(nz != NULL)
			*nz = _nz;
	}

	inline T& operator()(int x, int y, int z) {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
		if(!(z >= 0 && z < _nz)) throw "Index Z out of bounds!";
#endif 
		int index = lIndex(x,y,z);
		auto iter = siteMap.find(index);
		if(iter == siteMap.end())
			siteMap[index] = defaultValue;

		return siteMap[index];
	} 

	inline const T& operator()(int x, int y, int z) const {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
		if(!(z >= 0 && z < _nz)) throw "Index Z out of bounds!";
#endif 
		int index = lIndex(x,y,z);
		auto iter = siteMap.find(index);
		if(iter == siteMap.end())
			return defaultValue;
		return iter->second;
	}

	inline bool hasValue(const int &x, const int &y, const int &z){
		return siteMap.find(lIndex(x,y,z)) != siteMap.end();
	}

	inline int lIndex(const int &x, const int &y, const int &z) const {
		return z + _nz * (y + _ny*x);
	}
	/**
	 * Gets a pointer to the memory address
	 * that the actual array is located at.
	 * FFTW will need this eventually.
	 **/
	T* getArray() const {
		return _array;
	}

	int size() const {
		return _nx*_ny*_nz;
	}

};

}

#endif //_SPARSE_ARRAY3_H
