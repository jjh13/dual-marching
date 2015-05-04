/**
 * Generic class for handling 3-dimensional arrays, automatically
 * aligned for FFTW
 *
 * @author Joshua Horacsek
 *
 **/
#ifndef _ARRAY_H
#define _ARRAY_H

#include <string>
#include <memory>
#include <iostream>

// Useful for debugging
//#define NO_ARRAY_CHECK
namespace sisl {
template <class T, class Allocator = std::allocator<T> >
class array4 : private Allocator{

public:
	typedef T value_type;
	typedef typename Allocator::reference reference;
	typedef typename Allocator::const_reference const_reference;

	typedef typename Allocator::size_type size_type;
	typedef typename Allocator::difference_type difference_type;
	
	typedef typename Allocator::pointer iterator;
	typedef typename Allocator::const_pointer const_iterator;

	typedef Allocator allocator_type;
	allocator_type get_allocator() const {
		return static_cast<const Allocator&>(*this);
	}

	array4(const Allocator& a = Allocator()) : _nx(0), _ny(0), _nz(0), _nw(0), Allocator(a) {
		this->_array = NULL;
	}
	array4(const int &nx, const int &ny, const int &nz, const int &nw, const Allocator& a = Allocator()) :  _nx(nx), _ny(ny), _nz(nz), _nw(nw), Allocator(a) {
		int n = nx*ny*nz;
		size_type i;

		if(n == 0) {
			this->_array = nullptr;
			return;
		}
		this->first = Allocator::allocate(n);
		try {
			for (i = 0; i < n; ++i) {			
				Allocator::construct(first + i, T()); 
			}
		}catch(...) {
			for(size_type j = 0; j < i; ++j) {
		 		Allocator::destroy(first + j);
			}
			Allocator::deallocate(first, n);
			throw;
		}
		last = first + i;
		this->_array = static_cast<T *>(first);
	}

	~array4() {
		if (first != last && (_nx * _ny *_nz * _nw) != 0) {
			for (iterator i = first; i < last; ++i) {
				Allocator::destroy(i);
			}
			Allocator::deallocate(first, last - first);
		}
		this->_array = NULL;
	}



	inline array4& operator=(const array4 &rhs) {
		using namespace std; 

		int n = _nx*_ny*_nz;
		size_type i;

		if (this == &rhs)
			return *this;

		// Dealloc this object
		this->~array3();
		_nx = rhs._nx;
		_ny = rhs._ny;
		_nz = rhs._nz;
		_nw = rhs._nw;
		n = _nx*_ny*_nz*_nw;
		
		first = Allocator::allocate(n);
		if(n == 0) {
			this->_array = nullptr;
			return *this;
		}

		try {
			for (i = 0; i < n; ++i) {
				Allocator::construct(first + i, T(*(rhs.first+i))); 
			}
		}catch(...) {
			for(size_type j = 0; j < i; ++j) {
		 		Allocator::destroy(first + j);
			}
			Allocator::deallocate(first, n);
			throw;
		}
		last = first + i;
		this->_array = static_cast<T *>(first);
		return *this;
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

		return _array[w + _nw*(z + _nz * (y + _ny*x))];
	} 

	inline const T& operator()(int x, int y, int z, int w) const {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
		if(!(z >= 0 && z < _nz)) throw "Index Z out of bounds!";
		if(!(w >= 0 && w < _nw)) throw "Index W out of bounds!";
#endif 
		return _array[w + _nw*(z + _nz * (y + _ny*x))];
	}

	inline int lIndex(const int &x, const int &y, const int &z, const int &w) const {
		return w + _nw*(z + _nz * (y + _ny*x));
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

private:
	T *_array;
	unsigned int _nx, _ny, _nz, _nw;

	typename Allocator::pointer first;
	typename Allocator::pointer last;
};


template <class T, class Allocator = std::allocator<T> >
class array3 : private Allocator{

public:
	typedef T value_type;
	typedef typename Allocator::reference reference;
	typedef typename Allocator::const_reference const_reference;

	typedef typename Allocator::size_type size_type;
	typedef typename Allocator::difference_type difference_type;
	
	typedef typename Allocator::pointer       iterator;
	typedef typename Allocator::const_pointer const_iterator;

	typedef Allocator allocator_type;
	allocator_type get_allocator() const {
		return static_cast<const Allocator&>(*this);
	}

	array3(const Allocator& a = Allocator()) : _nx(0), _ny(0), _nz(0), Allocator(a) {
		this->_array = NULL;
	}
	array3(const int &nx, const int &ny, const int &nz, const Allocator& a = Allocator()) :  _nx(nx), _ny(ny), _nz(nz), Allocator(a) {
		int n = nx*ny*nz;

		size_type i;

		if(n == 0) {
			this->_array = nullptr;
			return;
		}
		this->first = Allocator::allocate(n);
		try {
			for (i = 0; i < n; ++i) {			
				Allocator::construct(first + i, T()); 
			}
		}catch(...) {
			for(size_type j = 0; j < i; ++j) {
		 		Allocator::destroy(first + j);
			}
			Allocator::deallocate(first, n);
			throw;
		}
		last = first + i;
		this->_array = static_cast<T *>(first);
	}

	~array3() {

		if (first != last && (_nx * _ny *_nz) != 0) {
			for (iterator i = first; i < last; ++i) {
				Allocator::destroy(i);
			}
			Allocator::deallocate(first, last - first);
		}
		this->_array = NULL;
	}



	inline array3& operator=(const array3 &rhs) {
		using namespace std; 

		int n = _nx*_ny*_nz;
		size_type i;

		if (this == &rhs)
			return *this;

		// Dealloc this object
		this->~array3();
		_nx = rhs._nx;
		_ny = rhs._ny;
		_nz = rhs._nz;
		n = _nx*_ny*_nz;
		
		first = Allocator::allocate(n);
		if(n == 0) {
			this->_array = nullptr;
			return *this;
		}

		try {
			for (i = 0; i < n; ++i) {
				Allocator::construct(first + i, T(*(rhs.first+i))); 
			}
		}catch(...) {
			for(size_type j = 0; j < i; ++j) {
		 		Allocator::destroy(first + j);
			}
			Allocator::deallocate(first, n);
			throw;
		}
		last = first + i;
		this->_array = static_cast<T *>(first);
		return *this;
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

		return _array[z + _nz * (y + _ny*x)];
	} 

	inline const T& operator()(int x, int y, int z) const {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
		if(!(z >= 0 && z < _nz)) throw "Index Z out of bounds!";
#endif 
		return _array[z + _nz * (y + _ny*x)];
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

private:
	T *_array;
	unsigned int _nx, _ny, _nz;

	typename Allocator::pointer first;
	typename Allocator::pointer last;
};

template <class T, class Allocator = std::allocator<T> >
class array2 : private Allocator{

public:
	typedef T value_type;
	typedef typename Allocator::reference reference;
	typedef typename Allocator::const_reference const_reference;

	typedef typename Allocator::size_type size_type;
	typedef typename Allocator::difference_type difference_type;
	
	typedef typename Allocator::pointer       iterator;
	typedef typename Allocator::const_pointer const_iterator;

	typedef Allocator allocator_type;
	allocator_type get_allocator() const {
		return static_cast<const Allocator&>(*this);
	}

	array2(const Allocator& a = Allocator()) : _nx(0), _ny(0), Allocator(a) {
		this->_array = NULL;
	}
	array2(const int &nx, const int &ny, const Allocator& a = Allocator()) :  _nx(nx), _ny(ny), Allocator(a) {
		int n = nx*ny;
		size_type i;

		if(n == 0) {
			this->_array = nullptr;
			return;
		}
		this->first = Allocator::allocate(n);
		try {
			for (i = 0; i < n; ++i) {			
				Allocator::construct(first + i, T()); 
			}
		}catch(...) {
			for(size_type j = 0; j < i; ++j) {
		 		Allocator::destroy(first + j);
			}
			Allocator::deallocate(first, n);
			throw;
		}
		last = first + i;
		this->_array = static_cast<T *>(first);
	}

	~array2() {

		if (first != last && (_nx * _ny) != 0) {
			for (iterator i = first; i < last; ++i) {
				Allocator::destroy(i);
			}
			Allocator::deallocate(first, last - first);
		}
		this->_array = NULL;
	}


	inline array2& operator=(const array2 &rhs) {
		
		int n = _nx*_ny;
		size_type i;

		if (this == &rhs)
			return *this;

		// Dealloc this object
		this->~array2();
		_nx = rhs._nx;
		_ny = rhs._ny;

		first = Allocator::allocate(n);

		if(n == 0) {
			this->_array = nullptr;
			return *this;
		}

		try {
			for (i = 0; i < n; ++i) {
				Allocator::construct(first + i, T(*(rhs.first+i))); 
			}
		}catch(...) {
			for(size_type j = 0; j < i; ++j) {
		 		Allocator::destroy(first + j);
			}
			Allocator::deallocate(first, n);
			throw;
		}
		this->_array = static_cast<T *>(first);
		return *this;
	}

	inline void getDims(unsigned int *nx, unsigned int *ny) const {
		if(nx != NULL)
			*nx = _nx;
		if(ny != NULL)
			*ny = _ny;
	}

	inline T& operator()(int x, int y) {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
#endif 

		return _array[(y + _ny*x)];
	} 

	inline const T& operator()(int x, int y) const {
#ifndef	NO_ARRAY_CHECK
		if(!(x >= 0 && x < _nx)) throw "Index X out of bounds!";
		if(!(y >= 0 && y < _ny)) throw "Index Y out of bounds!";
#endif 
		return _array[(y + _ny*x)];
	}

	inline int lIndex(const int &x, const int &y) const {
		return (y + _ny*x);
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
		return _nx*_ny;
	}

private:
	T *_array;
	unsigned int _nx, _ny;

	typename Allocator::pointer first;
	typename Allocator::pointer last;
};
}
#endif //_ARRAY_H
