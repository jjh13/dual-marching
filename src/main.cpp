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
#include <sisl/utility/scattered.hpp>
#include <sisl/utility/ply_writer.hpp>

#include <tuple>

#define VESION_STRING "0.1"

using namespace sisl;
using namespace std;
using namespace TCLAP;

template <class T>
class HamFunction {
public:
	HamFunction(){}	
	double evaluate(const double &x, const double &y, const double &z) {
		double xx = (x-0.5), yy = (y-0.5), zz = (z-0.5);
		return xx*xx + yy*yy + zz*zz - 0.25*0.25;

		// if( (x >= 0.25 && x <= 0.75) && (y >= 0.25 && y <= 0.75) && (z >= 0.25 && z <= 0.75)) {
		// 	return 1;
		// }
		// return -1.;
		//return -(sin(0.3141592654e1 * x) * sin(0.3141592654e1 * y) * sin(0.3141592654e1 * z) * (sqrt(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1)) - 0.2e1 * cos(0.8e1 * 0.3141592654e1 * (0.9e1 * z - 0.45e1) * pow(0.25e0 + pow(0.9e1 * x - 0.45e1, 0.2e1) + pow(0.9e1 * y - 0.45e1, 0.2e1) + pow(0.9e1 * z - 0.45e1, 0.2e1), -0.1e1 / 0.2e1)) - 0.2e1));
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

		ValueArg<std::string> outputArg("o", "output", "Output mesh/implicit file",true,"out.ply","filename");

		//cmd.add(outputArg);
		cmd.parse(argc, argv);
//		std::string output = outputArg.getValue();

		typedef linear_bcc_box<double, double> GFType;
		typedef bcc_odd<linear_bcc_box<double, double>, double, double> LATType;

		HamFunction<double> f;
	
		LATType *lat = new LATType(1./double(2.*100));

		lat->forEachLatticeSite([&](const int &i, const int &j, const int &k) {
			vector3<double> p = lat->getSitePosition(i,j,k);
			return f.evaluate(p);
		});

		dualConour<double>(lat);


		mc.marchLattice<LATType, double, double>(
			lat, 
			NULL, 
			NULL, 
			NULL, 
			0., 
			0.01, 
			vector3<double>(0.1,0.1,0.1), 
			vector3<double>(0.9,0.9,0.9));

		mc.writeSurface("marchingCubesTest.ply");

	}catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
	}catch (char const* e) {
		cerr << e << endl; 
	}
}

template <class T>
void dualConour(bcc_odd<linear_bcc_box<T,T>, T, T>  *lattice){
	int res = lattice->getResolution();
	T dh = lattice->getScale();


	std::vector<vector3<int> > vertices = {
		{ 0, 0, 0},
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

	std::vector<vector3<int> > centerIndices = {
		{-1, 0, 2}, 
		{0, -1, 2}, 
		{0, 1, 2}, 
		{1, 0, 2}, 
		{2, 0, 1}, 
		{2, 1, 0}, 
		{2,-1, 0}, 
		{2, 0, -1}, 
		{-2, 0, 1}, 
		{-2, 1, 0}, 
		{-2, -1, 0}, 
		{-2, 0, -1},
		{0, -2, 1}, 
		{1, -2, 0}, 
		{-1, -2, 0}, 
		{0, -2, -1}, 
		{0, 2, 1}, 
		{1, 2, 0},
		{-1, 2, 0}, 
		{0, 2, -1}, 
		{0, 1, -2}, 
		{-1, 0, -2}, 
		{1, 0, -2}, 
		{0, -1,-2}
	};

	std::vector<std::vector<int> > polygon_lookup = {
		{4 ,5  ,7  ,6},
		{8 ,9  ,10 ,11},
		{17,16 ,18 ,19},
		{12,13 ,15 ,14},
		{0 , 1 , 3 , 2},
		{20, 21, 22, 23},
	    
		{2 ,3 ,4 ,5 ,17,16},
		{5 ,17,19,20,22,7 },
		{1 ,3 ,4 ,6 ,13,12},
		{13,6 ,7 ,22,23,15},
		{0 ,2 ,16,18,9 ,8 },
		{18,9 ,11,21,20,19},
		{0 ,1 ,12,14,10,8 },
	    {10,14,15,23,21,11}
	};
	struct cell_vertex
	{
		cell_vertex(){

		}
		cell_vertex(T x, T y, T z) { 
			touching = {{x,y,z}};
		};
		std::vector<vector3<T> > touching;
		vector3<T> vertex;
		int vertexId;
	};

	std::vector<int> minimal_set = {1, 3, 5, 7, 8, 9, 11};
	sparse_array3<cell_vertex> face_hash_table(1000,1000,1000, {});
	std::vector<std::vector<vector3<int> > > faceList;

	utility::ply_writer<T> output_mesh;
	#pragma omp parallel for
	for(int i = 2; i < res*2-2; i+=2)
		for(int j = 2; j < res*2-2; j+=2)
			for(int k = 2; k < res*2-2; k++){
				int ii = i + (k%2);
				int jj = j + (k%2);
				int kk = k;

				// Get the value of this point
				T value = lattice->GV(ii, jj, kk);
				for(auto idx : minimal_set) {

					auto vdx = vertices[idx];
					int x = vdx.i + ii, 
						y = vdx.j + jj, 
						z = vdx.k + kk;
					T coeff = 0.5;
					T next_value = lattice->GV(x, y, z);
					vector3<T> pv;

					if((next_value > 0 && value > 0) || (next_value <0 && value < 0))
						continue;

					// Find the sign change.
					coeff = ((value - 0)/(value - next_value));
					pv = (vector3<T>(x,y,z) - vector3<T>(ii,jj,kk)).normalize();
					pv = vector3<T>(x,y,z) + pv * coeff;

					std::vector<int> luf = polygon_lookup[idx - 1];
					std::vector<vector3<int> > face;

					for(auto jdx : luf) {
						vector3<int> hash = centerIndices[jdx] + vector3<int>(ii*2, jj*2, kk*2);
						face.push_back(hash);

						{
							#pragma omp critical
							face_hash_table(hash.i, hash.j, hash.k).touching.push_back({hash.i, hash.j, hash.k});
						}
					}

					// Add the face to the list of faces
					{
						#pragma omp critical
						faceList.push_back(face);
					}
				}
			}

	// Calculate the vertex for each cell
	for (auto it = face_hash_table.siteMap.begin(); it != face_hash_table.siteMap.end(); ++it) {
		auto hash = it->first;
		auto vcache = it->second;
		auto pavg = vector3<T>(0,0,0);
		auto vList = vcache.touching;
		auto size = vList.size();


		for(auto v : vList) {
			pavg += v;
		}
		pavg = pavg * (1./((T)size))*dh;
		auto normal = lattice->grad_f(pavg).normalize();

		vertex3<T> vtx(pavg.i, pavg.j, pavg.k, normal.i, normal.j, normal.k);

		{
			face_hash_table.siteMap[hash].vertexId = output_mesh.addVertex(vtx);
		}
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


