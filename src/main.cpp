#include <iostream>
#include <tclap/CmdLine.h>

// SISL Main Include
#include <sisl/sisl.hpp>

// odd cartesian function spaces
#include <sisl/lattice/cartesian_odd.hpp>
#include <sisl/basis/tp3cubic.hpp>

// odd BCC function spaces
#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/basis/quintic.hpp>
#include <sisl/basis/linear_rdod.hpp>

// utility functions
#include <poisson/pointset.hpp>
#include <sisl/utility/isosurface.hpp>
#include <sisl/utility/dualbcc.hpp>
#include <sisl/utility/scattered.hpp>
#include <sisl/utility/ply_writer.hpp>

#include <Eigen/Dense>

#include <tuple>
#include <cmath>

#define VESION_STRING "0.1"

using namespace sisl;
using namespace std;
using namespace TCLAP;

template <class T>
class MarschnerLobb {
public:
	double a, f;

	MarschnerLobb(double Fm, double alpha) : f(Fm), a(alpha){ }
	double rho(double r){
		return cos(2*M_PI*f*cos(r*M_PI/2.));
	}	
	double evaluate(const double &xx, const double &yy, const double &zz) {
		double x = (xx-0.5)*2., y = (yy-0.5)*2., z = (zz-0.5)*2.;

		double r = rho(sqrt(x*x + y*y));


		double ret = 1. - sin(M_PI*z/2.) + a * (1. + r);

		return ret/(2. + 2.*a);

	}

	double evaluate(const vector3<T> &p) {return this->evaluate(p.i, p.j, p.k);}
	double div(const vector3<T> &p) { return this->div(p.i, p.j, p.k); }
	double div(const double &x, const double &y, const double &z) {return 0;}
	double laplacian(const vector3<T> &p){return this->laplacian(p.i, p.j, p.k);}
	double laplacian(const double &x, const double &y, const double &z) {return 0;}
};

template <class T>
class HamFunction {
public:
	HamFunction(){}
	double evaluate(const double &x, const double &y, const double &z) {
		double xx = (x-0.5), yy = (y-0.5), zz = (z-0.5);
		return xx*xx + yy*yy + zz*zz - 0.25*0.25;

		// if( (x >= 0.25 && x <= 0.75) && (y >= 0.25 && y <= 0.75) && (z >= 0.25 && z <= 0.75)) {
		// 	return -1;
		// }
		// return 1.;
		// return -(sin(0.3141592654e1 * x) * sin(0.3141592654e1 * y) * sin(0.3141592654e1 * z) * (sqrt(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1)) - 0.2e1 * cos(0.8e1 * 0.3141592654e1 * (0.9e1 * z - 0.45e1) * pow(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1), -0.1e1 / 0.2e1)) - 0.2e1));
	}

	double evaluate(const vector3<T> &p) {
		return this->evaluate(p.i, p.j, p.k);
	}
	
	double div(const vector3<T> &p) { return this->div(p.i, p.j, p.k); }
	double div(const double &x, const double &y, const double &z) {
		return 0;
	}
	
	double laplacian(const vector3<T> &p){return this->laplacian(p.i, p.j, p.k);}
	double laplacian(const double &x, const double &y, const double &z) {
		return 0;
	}
};

template <class T>
void dualConour(bcc_odd<linear_bcc_box<T,T>, T, T>  *lattice);

int main(int argc, char *argv[])
{
	try {
		CmdLine cmd("Dual BCC marching", ' ', VESION_STRING);
		sisl::utility::marchingCubes<double> mc;
		sisl::utility::dualbcc_isosurface<double> dbcc;

		ValueArg<std::string> outputArg("o", "output", "Output mesh/implicit file",true,"out.ply","filename");

		//cmd.add(outputArg);
		cmd.parse(argc, argv);
//		std::string output = outputArg.getValue();

		typedef linear_bcc_box<double, double> GFType;
		typedef bcc_odd<linear_bcc_box<double, double>, double, double> LATType;

		HamFunction<double> f;
		// MarschnerLobb<double> f(6, 0.25);

		LATType *lat = new LATType(1./double(2.*25));

		lat->forEachLatticeSite([&](const int &i, const int &j, const int &k) {
			vector3<double> p = lat->getSitePosition(i,j,k);
			return f.evaluate(p);// - 0.5;
		});

		dualConour<double>(lat);


		// dbcc.contour<LATType, double, double>(
		// 	lat, 0.5, lat->getScale(), 
		// 	vector3<double>(0.1,0.1,0.1), 
		// 	vector3<double>(0.9,0.9,0.9)
		// );

		// mc.marchLattice<LATType, double, double>(
		// 	lat, 
		// 	NULL, 
		// 	NULL, 
		// 	NULL, 
		// 	0, 
		// 	0.005, 
		// 	vector3<double>(0.1,0.1,0.1), 
		// 	vector3<double>(0.9,0.9,0.9));
		dbcc.writeSurface("dualBccTest.ply");
		// mc.writeSurface("marchingCubesTest.ply");

	}catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
	}catch (char const* e) {
		cerr << e << endl; 
	}
}

template<class T> 
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

	if(optimize && points.size() == 3){
		JacobiSVD<EMatrix> svd(A, ComputeThinU | ComputeThinV);
		EMatrix U = svd.matrixU(), V = svd.matrixV();
		EMatrix SS = EMatrix::Zero(V.cols(), U.rows());

		auto lambda = svd.singularValues();
		for(int i = 0; i < lambda.rows(); i++){
			if(abs(lambda(i)) > 0.1)
				SS(i,i) = 1./lambda(i);
		}
	}
	return center;
}

template <class T>
void dualConour(bcc_odd<linear_bcc_box<T,T>, T, T>  *lattice){
	int res = lattice->getResolution();
	T dh = lattice->getScale();
	res = int(1./dh)+1;


	std::vector<vector3<int> > vertices = {
	    {0,0,0},
	    {-1,0,1},
	    {0,1,1},
	    {1,0,1},
	    {0,-1,1},
	    {-1,-1,0},
	    {1,-1,0},
	    {1,1,0},
	    {-1,1,0},
	    {0,-1,-1},
	    {1,0,-1},
	    {0,1,-1},
	    {-1,0,-1}
	};
	std::vector<vector3<int> > centerIndices = {
	    { 2, 0, 0},
	    {-2, 0, 0},
	    { 0, 2, 0},
	    { 0,-2, 0},
	    { 0, 0, 2},
	    { 0, 0,-2},
	    { 1, 1, 1},
	    { 1, 1,-1},
	    { 1,-1, 1},
	    { 1,-1,-1},
	    {-1, 1, 1},
	    {-1, 1,-1},
	    {-1,-1, 1},
	    {-1,-1,-1},
	};

	std::vector<std::vector<int>> adj_index = {

	    {4, 10, 12, 1},
	    {4, 10, 6, 2},
	    {4, 6, 8, 0},
	    {4, 8, 12, 3},
	    {3, 13, 12, 1},
	    {0, 8, 9, 3},
	    {2,6,7,0},
	    {1,10,11,2},
	    {3,9,13, 5},
	    {0,7,9,5},
	    {2,11,7, 5},
	    {1,13,11,5}
	};

	std::vector<std::vector<std::vector<int>>> polygon_lookup = {
	    {{4, 10, 12}, {10, 12, 1}},
	    {{4, 10, 6}, {10, 6, 2}},
	    {{4, 6, 8}, {6, 8, 0}},
	    {{4, 8, 12}, {8, 12, 3}},
	    {{3, 13, 12}, {8, 12, 3}},
	    {{0, 8, 9}, {8, 9, 3}},
	    {{2,6,7}, {6,7,0}},
	    {{1,10,11}, {10,11,2}},
	    {{3,9,13}, {9,13,5}},
	    {{0,7,9}, {7,9,5}},
	    {{2,11,7}, {11,7,5}},
	    {{1,13,11}, {13,11,5}}
	};

	struct cell_vertex
	{
		cell_vertex(){}
		std::vector<vector3<T> > touching;
		vector3<T> vertex;
		int vertexId;
	};

	std::vector<int> minimal_set = {1 , 2, 3, 4, 5 };
	sparse_array3<cell_vertex> face_hash_table(1000,1000,1000, {});
	std::vector<std::vector<vector3<int> > > faceList;

	utility::ply_writer<T> output_mesh;
	#pragma omp parallel for
	for(unsigned int i = 2; i < res-2; i+=2)
		for(unsigned int j = 2; j < res - 2; j++)
			for(unsigned int k = 2; k < res - 2; k++) {

				int ii = i + ((j&1) ^ (k&1));

				int jj = j;
				int kk = k;
				// if(j % 2) {
				// 		if(k % 2){
				// 			ii = i*2;
				// 		}else{
				// 			ii = i*2 + 1;
				// 		}
				// 	} else {
				// 		if(k % 2){
				// 			ii = i*2 + 1;
				// 		}else{
				// 			ii = i*2;
				// 		}
				// 	}

				// Get the value of this point
				T value = lattice->f(ii*dh, jj*dh, kk*dh);
				for(auto idx : minimal_set) {

					auto vdx = vertices[idx];
					int x = vdx.i + ii, 
						y = vdx.j + jj, 
						z = vdx.k + kk;
					T coeff = 0.5;
					T next_value = lattice->f(x*dh, y*dh, z*dh);
					vector3<T> pv, n;

					if((next_value > 0 && value > 0) || (next_value <0 && value < 0))
						continue;

					// Find the sign change.
					coeff = ((value - 0)/(value - next_value));
					pv = n = vector3<T>(x,y,z) - vector3<T>(ii,jj,kk);
					pv = vector3<T>(ii,jj,kk) + pv * coeff;
					n = n * (value > next_value ? -1 : 1);

					/*  */
					std::vector<int> adj = adj_index[idx - 1];
					std::vector<std::vector<int>> luf = polygon_lookup[idx - 1];

					// Build triangles
					for(auto triangle : luf) {
						std::vector<vector3<int> > tri;

						vector3<int> 
								hash1 = centerIndices[triangle[0]] + vector3<int>(ii*2, jj*2, kk*2),
								hash2 = centerIndices[triangle[1]] + vector3<int>(ii*2, jj*2, kk*2),
								hash3 = centerIndices[triangle[2]] + vector3<int>(ii*2, jj*2, kk*2);

						vector3<int> t = (hash2 - hash1)%(hash3 - hash1);
						vector3<T> dir(t.i, t.j, t.k);
						if(dir * n > 0) tri = {hash1, hash2, hash3};
						else tri = {hash3, hash2, hash1};

						// 
						{
							#pragma omp critical
							faceList.push_back(tri);
						}
					}
					
					// Mark the hashed vertex as having seen this vertex
					for(auto jdx : adj) {
						vector3<int> hash = centerIndices[jdx] + vector3<int>(ii*2, jj*2, kk*2);
						{
							#pragma omp critical
							face_hash_table(hash.i, hash.j, hash.k).touching.push_back({pv.i, pv.j, pv.k});
						}
					}
				}
			}

	// Calculate the vertex for each cell
	for (auto it = face_hash_table.siteMap.begin(); it != face_hash_table.siteMap.end(); ++it) {
		auto hash = it->first;
		auto vcache = it->second;

		std::vector<vector3<T>> normals;
		for(auto v : vcache.touching) 
			normals.push_back(lattice->grad_f(v*dh).normalize());
		

		auto pavg = optimize_for_feature<T>(vcache.touching, normals) * dh;
		auto normal = lattice->grad_f(pavg).normalize();
		
		face_hash_table.siteMap[hash].vertexId = output_mesh.addVertex({pavg, normal});
	}

	/* Build the final face list */
	for(auto face : faceList) {
		std::vector<int> index_face; 
		for(auto hash : face) {
			index_face.push_back(face_hash_table(hash.i, hash.j, hash.k).vertexId);
		}
		output_mesh.addPolygon(index_face);
	}

	output_mesh.writePly("output.ply");
}


