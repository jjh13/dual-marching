#ifndef _PLY_WRITER_H_
#define _PLY_WRITER_H_

#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>

#include <sisl/sisl.hpp>

namespace sisl{
namespace utility{

template <class T>
class ply_writer {
public:
	ply_writer() : currentIndex(0){}
	~ply_writer(){}

	bool writePly(const std::string &out) const {
		using namespace std;
		using namespace sisl;

		ofstream fp(out.c_str(), ios::out | ios::binary);
		if(!fp.good()) return false;

		fp << "ply" << endl;
		fp << "format binary_little_endian 1.0" << endl;
		fp << "element vertex " << verts.size() << endl;
		fp << "property float x" << endl;
		fp << "property float y" << endl;
		fp << "property float z" << endl;
		fp << "property float nx" << endl;
		fp << "property float ny" << endl;
		fp << "property float nz" << endl;
		fp << "element face " << faces.size()+polygons.size() << endl;
		fp << "property list uchar int vertex_indices" << endl;
		fp << "end_header" << endl;

		for(unsigned int i = 0; i < verts.size(); i++) {
			struct _p {
				float x, y, z, i, j , k;
			}__attribute__ ((packed)) p;

			p.x = verts[i].p.i;
			p.y = verts[i].p.j;
			p.z = verts[i].p.k;
			p.i = verts[i].n.i;
			p.j = verts[i].n.j;
			p.k = verts[i].n.k;

			//printf("v %f  %f  %f\nvn %f %f %f \n", p.x, p.y, p.z, p.i, p.j, p.k);
			fp.write((char*)&p, sizeof(_p));
		}

		for(unsigned int i = 0; i < faces.size(); i++) {
			struct _f {
				char t;
				int i,j,k;
			} __attribute__ ((packed)) f;

			f.t = 3;
			f.i = std::get<0>(faces[i]);
			f.j = std::get<1>(faces[i]);
			f.k = std::get<2>(faces[i]);
			fp.write((char*)&f, sizeof(_f));
		}

		for(unsigned int i = 0; i < polygons.size(); i++) {
			auto poly = polygons[i];

			char t = (char)poly.size();
			fp.write((char*)&t, sizeof(char));
			for(int index : poly)
				fp.write((char*)&index, sizeof(int));
		}

		fp.close();
		return true;
	}

	bool transformMesh(const sisl::matrix4x4<T> &t) {
		for(unsigned int i = 0; i < verts.size(); i++) 
			verts[i].p = (t * verts[i].p);
	}

	// Returns the index of the new vertex
	int addVertex(const sisl::vertex3<T> &p) {
		using namespace std;
		using namespace sisl;
		int r = 0;

		#pragma omp critical (add_vertex)
		{
			verts.push_back(p);
			r = (currentIndex++);
		}

		return r;
	}

	int addTriangle(sisl::vertex3<T> v1, sisl::vertex3<T> v2, sisl::vertex3<T> v3) {
		std::tuple<int,int,int> f;
		int r = 0;
		#pragma omp critical (add_face)
		{
			int i1 = addVertex(v1),
				i2 = addVertex(v2),
				i3 = addVertex(v3);
			
			r = addPolygon({i1,i2,i3});
		}
		return r; 
	}

	int addTriangle(int i, int j, int k) {
		return addPolygon({i,j,k});
	}
	
	int addPolygon(const std::vector<int> &v) {
		int r = 0;

		#pragma omp critical (add_face) 
		{
			polygons.push_back(v);
			r = int(polygons.size())-1;
		}
		return r;
	}

	int countVertices() {
		return currentIndex;
	}

	int countFaces() {
		return int(faces.size()+polygons.size());
	}

private:
	int currentIndex;

	// This maps verticies to their indices
	std::multimap<sisl::vertex3<T>, int> vertexIndexMap;
	std::vector<sisl::vertex3<T>> verts;
	std::vector<std::tuple<int,int,int>> faces;

	std::vector<std::vector<int>> polygons;
};
};
};

#endif // _PLY_OBJECT_H_