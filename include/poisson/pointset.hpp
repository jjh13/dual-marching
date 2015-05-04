#ifndef _POINT_SET_H_
#define _POINT_SET_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sisl/sisl.hpp>

template <class T>
class pointset {
public:
	bool readXyz(const std::string &in) {
		std::ifstream fp(in.c_str(), std::ios::in | std::ios::binary);
		std::string line;

		if(!fp.good()) return false;

		while(!fp.eof()) {
			std::string line;
			getline(fp, line);
			char *tokens = strtok((char*)line.c_str(), " \t");
			if(tokens == NULL)
				continue;
			T x = atof(tokens);
			T y = atof(strtok(NULL, " \t"));
			T z = atof(strtok(NULL, " \t"));
			T i = atof(strtok(NULL, " \t"));
			T j = atof(strtok(NULL, " \t"));
			T k = atof(strtok(NULL, " \t"));
			points.push_back(sisl::vertex3<T>(x,y,z,i,j,k));
		}
		fp.close();
		return true;
	}

	std::vector<sisl::vector3<T>> getPositionVector(){
		std::vector<sisl::vector3<T>> l;
		for(unsigned int i = 0; i < points.size(); i++) 
			l.push_back(points[i].p);
		return l;
	};

	std::vector<sisl::vector3<T>> getNormalVector(){
		std::vector<sisl::vector3<T>> l;
		for(unsigned int i = 0; i < points.size(); i++) 
			l.push_back(points[i].n);
		return l;
	};

	std::vector<T> getWeigthedNormalColumnVector(sisl::vector3<T> w){
		std::vector<T> l;
		for(unsigned int i = 0; i < points.size(); i++) 
			l.push_back(points[i].n * w);
		return l;
	};
	
	void applyTransform(const sisl::matrix4x4<T> &transform) {
		for(int i = 0; i < points.size(); i++)
			points[i].p = transform * points[i].p;
	}

	void applyNormalTransform(const sisl::matrix4x4<T> &transform) {
		for(int i = 0; i < points.size(); i++)
			points[i].n = transform * points[i].n;
	}

	sisl::vector3<T> getCenter() {
		sisl::vector3<T> c;
		for(int i = 0; i < points.size(); i++)
			c += points[i].p;	
		return c*(1./double(points.size()));
	}

	T calculateRadius(sisl::vector3<T> center) {
		T scale = 0.;
		for(int i = 0; i < points.size(); i++) {
			T distance = (center - points[i].p).norm();
			if(distance > scale)
				scale = distance;
		}
		return scale;
	}

	int count() {
		return points.size();
	} 
private:
	std::vector<sisl::vertex3<T>> points;
};

#endif //_POINT_SET_H_