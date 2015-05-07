#ifndef _DUAL_BCC_ISOSURFACE_H_
#define _DUAL_BCC_ISOSURFACE_H_

#include <sisl/sisl.hpp>
#include <sisl/sparse_array.hpp>
#include <sisl/utility/ply_writer.hpp>

#include <Eigen/Dense>

#include <tuple>
#include <unordered_map>
#include <map>

namespace sisl{
namespace utility{
using namespace std;

template<class T>
class dualbcc_isosurface{
public:
	struct cell_vertex
	{
		cell_vertex(){}
		std::vector<vector3<T> > touching;
		vector3<T> vertex;
		int vertexId;
	};

	dualbcc_isosurface() : face_hash_table(1000,1000,1000, {}){
		this->faceList.clear();
	}

	template<class L, class I, class O>
	void contour(
			L *f,
			const O &isoValue, 
			const I &scalingParameter,
			sisl::vector3<I> origin,
			sisl::vector3<I> boundary ){

		I dh = scalingParameter;
		int res = int(1./(2.*dh));

		// Go over every lattice point
		#pragma omp parallel for
		for(int i = 2; i < res*2-2; i+=2) {
			std::vector<std::vector<vector3<int>>> localFaceList;

			for(int j = 2; j < res*2-2; j+=2)
				for(int k = 2; k < res*2-2; k++){
					int ii = i + (k%2);
					int jj = j + (k%2);
					int kk = k;

					O value = f->f(dh*ii, dh*jj, dh*kk) - isoValue;

					// For each face in the minimal amount of faces of
					// the polyhedron
					for(auto idx : minimal_face_set) {
						auto polyhedron_vertex = polyhedron_vertices[idx];
						int x = polyhedron_vertex.i + ii, 
							y = polyhedron_vertex.j + jj, 
							z = polyhedron_vertex.k + kk;
						O next_value = f->f(dh*x, dh*y, dh*z) - isoValue;
						I zero_solution = 0.5; 

						vector3<T> pv, n;

						// No sign change?
						if((next_value > 0 && value > 0) || (next_value <0 && value < 0))
							continue; // Whatever

						// Find the sign change.
						zero_solution = ((value - 0)/(value - next_value));
						pv = n = vector3<T>(x,y,z) - vector3<T>(ii,jj,kk);
						pv = vector3<T>(ii,jj,kk) + pv * zero_solution;
						n = n * (value > next_value ? -1 : 1);


						// Lookup all the dual points that touch this vertex
						std::vector<int> adj = adj_index[idx - 1];
						std::vector<std::vector<int>> luf = triangle_lookup[idx - 1];

						// Push all the faces into our local face list.
						for(auto triangle : luf) {
							vector3<int> 
									hash1 = center_hash_offsets[triangle[0]] + vector3<int>(ii*2, jj*2, kk*2),
									hash2 = center_hash_offsets[triangle[1]] + vector3<int>(ii*2, jj*2, kk*2),
									hash3 = center_hash_offsets[triangle[2]] + vector3<int>(ii*2, jj*2, kk*2);

							vector3<int> t = (hash2 - hash1)%(hash3 - hash1);
							vector3<T> dir(t.i, t.j, t.k);
							if(dir * n > 0) localFaceList.push_back((std::vector<vector3<int>>){hash1, hash2, hash3});
							else localFaceList.push_back((std::vector<vector3<int>>){hash3, hash2, hash1});
						}
						
						// Mark the hashed dual vertex as having seen this primal vertex
						for(auto jdx : adj) {
							vector3<int> hash = center_hash_offsets[jdx] + vector3<int>(ii*2, jj*2, kk*2);
							
							#pragma omp critical (hash_bash_bcc)
							{
								face_hash_table(hash.i, hash.j, hash.k).touching.push_back({pv.i, pv.j, pv.k});
							}
						}

					}
				}
				// Merge the faces back in to global face list
				#pragma omp critical (lizst_cyst)
				{
					faceList.reserve(faceList.size() + localFaceList.size());
					faceList.insert(faceList.end(), localFaceList.begin(), localFaceList.end());
				}
		}
		processVertices<L,I,O>(f, dh);
		processFaces();
	}

	bool writeSurface(const std::string &out) const {
		return output_mesh.writePly(out);
	}

private:
	template<class L, class I, class O>
	void processVertices(L *f, const I &dh){
		// Calculate the vertex for each cell
		for (auto it = face_hash_table.siteMap.begin(); it != face_hash_table.siteMap.end(); ++it) {
			auto hash = it->first;
			auto vcache = it->second;

			std::vector<vector3<T>> normals;
			for(auto v : vcache.touching) 
				normals.push_back(f->grad_f(v*dh).normalize());
			

			auto pavg = optimize_for_feature(vcache.touching, normals) * dh;
			auto normal = f->grad_f(pavg).normalize();
			
			face_hash_table.siteMap[hash].vertexId = output_mesh.addVertex({pavg, normal});
		}
	}

	void processFaces(){
		/* Build the final face list */
		for(auto face : faceList) {
			std::vector<int> index_face; 
			for(auto hash : face) {
				index_face.push_back(face_hash_table(hash.i, hash.j, hash.k).vertexId);
			}
			output_mesh.addPolygon(index_face);
		}
	}

	std::vector<std::vector<vector3<int>>> faceList;
	sisl::sparse_array3<cell_vertex> face_hash_table; 
	utility::ply_writer<T> output_mesh;
	
	const std::vector<vector3<int>> polyhedron_vertices = {
		{ 0, 0, 0},	{ 2, 0, 0}, {-2, 0, 0}, { 0, 2, 0},	{ 0,-2, 0},
		{ 0, 0, 2},	{ 0, 0,-2}, { 1, 1, 1},	{ 1, 1,-1},	{ 1,-1, 1},	
		{ 1,-1,-1}, {-1, 1, 1},	{-1, 1,-1},	{-1,-1, 1},	{-1,-1,-1},
	};

	const std::vector<vector3<int>> center_hash_offsets = {
		{-1, 0, 2}, {0, -1, 2}, {0, 1, 2}, 	{1, 0, 2}, 
		{2, 0, 1}, 	{2, 1, 0}, {2,-1, 0}, 	{2, 0, -1}, 
		{-2, 0, 1},	{-2, 1, 0},	{-2, -1, 0}, {-2, 0, -1},
		{0, -2, 1}, {1, -2, 0}, {-1, -2, 0}, {0, -2, -1}, 
		{0, 2, 1}, 	{1, 2, 0}, {-1, 2, 0},	{0, 2, -1}, 
		{0, 1, -2},	{-1, 0, -2},  {1, 0, -2}, {0, -1,-2}
	};

	const std::vector<std::vector<int>> adj_index = {
		{4 ,5  ,7  ,6}, {8 ,9  ,10 ,11},
		{17,16 ,18 ,19}, {12,13 ,15 ,14},
		{0 , 1 , 3 , 2}, {20, 21, 22, 23},
		{2 ,3 ,4 ,5 ,17,16}, {5 ,17,19,20,22,7 },
		{1 ,3 ,4 ,6 ,13,12}, {13,6 ,7 ,22,23,15},
		{0 ,2 ,16,18,9 ,8 }, {18,9 ,11,21,20,19},
		{0 ,1 ,12,14,10,8 }, {10,14,15,23,21,11}
	};

	const std::vector<std::vector<std::vector<int>>> triangle_lookup = {
		{{4,5,7}, {4,7,6}},	{{8,9,10}, {8, 10, 11}},	
		{{17,16,18}, {17, 18, 19}}, {{12, 13, 15}, {12, 15, 14}},
		{{0,1,3}, {0,3,2}},	{{20,21,22}, {20,22,23}}, 
		{{2,3,4}, {2,4,5}, {2,5,17}, {2,17,16}}, {{5,17,19}, {5,19,20}, {5,20,22}, {5, 22, 7}}, 
		{{1,3,4}, {1,4,6}, {1,6,13}, {1,13,12}}, {{13, 6, 7}, {13, 7,22}, {13,22, 23}, {13, 23, 15}}, 
		{{0,2,16}, {0, 16, 18}, {0, 18,9}, {0,9,8}}, {{18,9,11}, {18, 11, 21}, {18, 21,20}, {18,20,19}}, 
		{{0,1,12}, {0,12,14}, {0, 14, 10}, {0, 10, 8}}, {{10,14,15}, {10, 15, 23}, {19,23,21}, {10,21,11}} 
	};

	const std::vector<int> minimal_face_set = {1, 3, 5, 7, 8, 9, 11};


	vector3<T> optimize_for_feature(
			const std::vector<vector3<T>> &points, 
			const std::vector<vector3<T>> &normals,
			const T &threshold = 0.1,
			const bool &optimize = true) {
		using namespace Eigen;
		typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> EMatrix;
		vector3<T> center(0,0,0);
		EMatrix A(points.size(), 3), b(points.size(), 1);
		unsigned int i = 0; 

		// Calculate the center and setup the matix
		for(auto v : points) { 
			center += v; 

			A(i, 0) = v.i;
			A(i, 1) = v.j;
			A(i, 2) = v.k;

			b(i, 0) = points[i] * normals[i];
			i++;
		}

		center = center * (1./(T(points.size())));

		// if(optimize && points.size() == 3){
		// 	JacobiSVD<EMatrix> svd(A, ComputeThinU | ComputeThinV);
		// 	EMatrix U = svd.matrixU(), V = svd.matrixV();
		// 	EMatrix SS = EMatrix::Zero(V.cols(), U.rows());

		// 	auto lambda = svd.singularValues();
		// 	for(int i = 0; i < lambda.rows(); i++){
		// 		if(abs(lambda(i)) > 0.1)
		// 			SS(i,i) = 1./lambda(i);
		// 	}

		// 	cout << "A" << endl;
		// 	cout << A << endl;

		// 	cout << "U" << endl;
		// 	cout << U << endl;

		// 	cout << "SS" << endl;
		// 	cout << lambda << endl;

		// 	cout << "V" << endl;
		// 	cout << V << endl;

		// 	cout << "LS" << endl;
		// 	cout << SS << endl;

		// 	cout << "B" << endl;
		// 	EMatrix sol = (V.transpose()*SS.transpose()*U.transpose() )*b;// << endl;
		// 	cout << sol << endl;

		// 	sol = svd.solve(b);
		// 	cout << "B" << endl;
		// 	cout << sol << endl;

		// //	center = vector3<T>(sol(0,0), sol(1,0), sol(2,0));


		// }


		return center;
	}

};
};
};

#endif // _DUAL_BCC_ISOSURFACE_H_