#define VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENMP
#define VIENNACL_HAVE_EIGEN

#include <typeinfo>

#include <sisl/sisl.hpp>
#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/lattice/cartesian_odd.hpp>
#include <sisl/lattice/fcc_odd.hpp>

#include <Eigen/Sparse>

#include "viennacl/scalar.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"


namespace sisl{
namespace utility{
extern const char * _bcc_optimization_kernel_dbl;
extern const char * _cc_optimization_kernel_dbl;

/**
 * Samples the basis function specified by 
 * GF at locations specified by samplePts,
 * returns an Eigen::SparseMatrix that 
 * corresponds to \mathbf{Phi} in the paper
 */
template <class L, class GF, class I, class O>
void createBigPhi(
		const std::vector<vector3<I>> &samplePts, 
		L *l,
		vector3<I> latticeShift,
		viennacl::compressed_matrix<O> &phi,
		viennacl::compressed_matrix<O> &phiT) {

	int nSize = l->numberOfLatticeSites();
	std::vector<Eigen::Triplet<O>> eigenMap;
	std::vector<std::tuple<int,int,int>> st = GF::getSupport();
	I dh = l->getScale();

	// Sample the basis function to create the Xi vector
	for(unsigned int ptIndex = 0; ptIndex < samplePts.size(); ptIndex++) {
		vector3<I> pt = samplePts[ptIndex] + (latticeShift * dh);
		vector3<int> pIdx = l->getNearestIndex(pt);

		printf("%f %f %f;\r",pt.i, pt.j,pt.k );
		// Sample each lattice site
		for(unsigned int stIdx = 0; stIdx < st.size(); stIdx++) {
			std::tuple<int,int,int> s = st[stIdx];
			vector3<I> q = pt - l->getSitePosition(
				std::get<0>(s) + pIdx.i,
				std::get<1>(s) + pIdx.j,
				std::get<2>(s) + pIdx.k);

			int mindex = l->lIndex(
				std::get<0>(s) + pIdx.i,
				std::get<1>(s) + pIdx.j,
				std::get<2>(s) + pIdx.k);
			
			O value = GF::M(q.i/dh, q.j/dh, q.k/dh);

			if(value == 0 || mindex < 0) continue;
			eigenMap.push_back(Eigen::Triplet<O>(ptIndex, mindex, value));
		}
	}
	Eigen::SparseMatrix<O> cpuPhi(samplePts.size(), nSize);
	Eigen::SparseMatrix<O> cpuPhiT(nSize, samplePts.size());

	cpuPhi.setFromTriplets(eigenMap.begin(), eigenMap.end());
	cpuPhiT = Eigen::SparseMatrix<O>(cpuPhi.transpose());

	printf("\nCopying to gpu...\n");
	viennacl::copy(cpuPhi, phi);
	viennacl::copy(cpuPhiT, phiT);

	printf("Copying OK...\n");
}


template <class L, class GF, class I, class O>
void splatSamplesOntoLattice(
		const std::vector<vector3<I>> &samplePts,
		const std::vector<O> &samples,
		L *l,
		vector3<I> latticeShift) {

	int nSize = l->numberOfLatticeSites();
	std::vector<Eigen::Triplet<O>> eigenMap;
	std::vector<std::tuple<int,int,int>> st = GF::getSupport();
	I dh = l->getScale();

	// Sample the basis function to create the Xi vector
	for(unsigned int ptIndex = 0; ptIndex < samplePts.size(); ptIndex++) {
		vector3<I> pt = samplePts[ptIndex] + (latticeShift * dh);
		vector3<int> pIdx = l->getNearestIndex(pt);

		// Sample each lattice site
		//#pragma omp parallel for
		for(unsigned int stIdx = 0; stIdx < st.size(); stIdx++) {
			std::tuple<int,int,int> s = st[stIdx];
			vector3<I> q = pt - l->getSitePosition(
				std::get<0>(s) + pIdx.i,
				std::get<1>(s) + pIdx.j,
				std::get<2>(s) + pIdx.k);

			O value = GF::M(q.i/dh, q.j/dh, q.k/dh) * samples[ptIndex] + l->GV(
					std::get<0>(s) + pIdx.i,
					std::get<1>(s) + pIdx.j,
					std::get<2>(s) + pIdx.k);

			l->SV(
				std::get<0>(s) + pIdx.i,
				std::get<1>(s) + pIdx.j,
				std::get<2>(s) + pIdx.k, value);
		}
	}
}

template <class L, class GF, class I, class O>
bool ConjugateGradientBCC(
		viennacl::compressed_matrix<O> &phi,
		viennacl::compressed_matrix<O> &phiT,
		const std::vector<O> &samples,
		const std::vector<std::tuple<int,int,int,O>> &optimizationVector,
		L *l,
		O *err = NULL, int *iter = NULL) {


	int nSize = l->numberOfLatticeSites();
	int niter = 0;

/*	std::vector<double> cpuX(nSize);
	l->forEachLatticeSite([&](const int &i, const int &j, const int &k) {
		return cpuX[l->lIndex(i,j,k)] = l->GV(i,j,k);
	});*/

	// Setup the residual vectors
	viennacl::vector<O> residual(nSize), p(nSize), Ap(nSize), x(nSize);
	viennacl::vector<O> apTT(samples.size()), lsample(samples.size());


	// Initialize with zero 
	residual.clear(); p.clear();
	Ap.clear(); x.clear(); apTT.clear();

//	viennacl::copy(cpuX, x);

	O resNew = 0;
	O resOld = 0;
	O alpha = 0;

	viennacl::copy(samples, lsample);

	// Do an initial iteration
	residual = viennacl::linalg::prod(phiT, lsample);
	p = viennacl::linalg::prod(phiT, lsample);
	resNew = viennacl::linalg::inner_prod(p,p);
	resOld = resNew;

	viennacl::ocl::handle<cl_mem> xr, yr, zr, blw;
	viennacl::ocl::program &normMtxPrgm = viennacl::ocl::current_context().add_program(sisl::utility::_bcc_optimization_kernel_dbl, "blCode");

	// Setup the opencl kernel for execution
	{
		std::vector<cl_int> xv, yv, zv;
		std::vector<cl_double> w;

		for(auto itr = optimizationVector.begin(); itr != optimizationVector.end(); ++itr) {
			auto t = (*itr);
			xv.push_back(std::get<0>(t));
			yv.push_back(std::get<1>(t));
			zv.push_back(std::get<2>(t));
			w.push_back(std::get<3>(t));
		}

		xr = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_int), &(xv[0]));

		yr = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_int), &(yv[0]));

		zr = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_int), &(zv[0]));

		blw = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_double), &(w[0]));
	}

	viennacl::ocl::kernel &blKern = normMtxPrgm.get_kernel("blnormbcc");
	//normMtxPrgm.add_kernel(blKern, "blnormbcc");

	for(niter = 0; niter < 5000; niter++) {
		alpha = 0;
		apTT = viennacl::linalg::prod(phi, p);
		Ap = viennacl::linalg::prod(phiT, apTT);

		viennacl::ocl::enqueue(
			blKern(p, 
				xr, yr, zr, blw, Ap,
				static_cast<cl_uint>(nSize),
				static_cast<cl_uint>(l->getResolution()),
				static_cast<cl_double>(1)));
		viennacl::backend::finish();

		alpha = resOld/viennacl::linalg::inner_prod(p, Ap);
		x = x + alpha * p;
		residual = residual - (alpha * Ap);
		resNew = viennacl::linalg::inner_prod(residual, residual);

		if(sqrt(resNew) < 1e-5)
			break;

		printf("%d: %f...\r",niter, sqrt(resNew));
		std::cout.flush();

		p = residual + (resNew/resOld)*p;
		resOld=resNew;
	}

	residual.clear();
	p.clear(); 
	Ap.clear();
	apTT.clear();
	
	std::vector<double> cpuX(nSize);
	viennacl::copy(x, cpuX);
	x.clear();

	l->forEachLatticeSite([&](const int &i, const int &j, const int &k) {
		return cpuX[l->lIndex(i,j,k)];
	});
	return true;
}

template <class L, class GF, class I, class O>
bool ConjugateGradientCC(
		viennacl::compressed_matrix<O> &phi,
		viennacl::compressed_matrix<O> &phiT,
		const std::vector<O> &samples,
		const std::vector<std::tuple<int,int,int,O>> &optimizationVector,
		L *l,
		O *err = NULL, int *iter = NULL) {


	int nSize = l->numberOfLatticeSites();
	int niter = 0;

	// Setup the residual vectors
	viennacl::vector<O> residual(nSize), p(nSize), Ap(nSize), x(nSize);
	viennacl::vector<O> apTT(samples.size()), lsample(samples.size());


	// Initialize with zero 
	residual.clear(); p.clear();
	Ap.clear(); x.clear(); apTT.clear();

	O resNew = 0;
	O resOld = 0;
	O alpha = 0;

	viennacl::copy(samples, lsample);

	// Do an initial iteration
	residual = viennacl::linalg::prod(phiT, lsample);
	p = viennacl::linalg::prod(phiT, lsample);
	resNew = viennacl::linalg::inner_prod(p,p);
	resOld = resNew;

	viennacl::ocl::handle<cl_mem> xr, yr, zr, blw;
	viennacl::ocl::program &normMtxPrgm = viennacl::ocl::current_context().add_program(sisl::utility::_cc_optimization_kernel_dbl, "blCode");

	// Setup the opencl kernel for execution
	{
		std::vector<cl_int> xv, yv, zv;
		std::vector<cl_double> w;

		for(auto itr = optimizationVector.begin(); itr != optimizationVector.end(); ++itr) {
			auto t = (*itr);
			xv.push_back(std::get<0>(t));
			yv.push_back(std::get<1>(t));
			zv.push_back(std::get<2>(t));
			w.push_back(std::get<3>(t));
		}

		xr = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_int), &(xv[0]));

		yr = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_int), &(yv[0]));

		zr = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_int), &(zv[0]));

		blw = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						1024 * sizeof(cl_double), &(w[0]));
	}

	viennacl::ocl::kernel &blKern = normMtxPrgm.get_kernel("blnorm");
	//normMtxPrgm.add_kernel(blKern, "blnorm");

	for(niter = 0; niter < 5000; niter++) {
		alpha = 0;
		apTT = viennacl::linalg::prod(phi, p);
		Ap = viennacl::linalg::prod(phiT, apTT);

		viennacl::ocl::enqueue(
			blKern(p, 
				xr, yr, zr, blw, Ap,
				static_cast<cl_uint>(nSize),
				static_cast<cl_uint>(l->getResolution()),
				static_cast<cl_uint>(l->getResolution()),
				static_cast<cl_uint>(l->getResolution()),
				static_cast<cl_double>(1)));
		viennacl::backend::finish();

		alpha = resOld/viennacl::linalg::inner_prod(p, Ap);
		x = x + alpha * p;
		residual = residual - (alpha * Ap);
		resNew = viennacl::linalg::inner_prod(residual, residual);

		if(sqrt(resNew) < 1e-5)
			break;

		printf("%d: %f...\r",niter, sqrt(resNew));
		std::cout.flush();

		p = residual + (resNew/resOld)*p;
		resOld=resNew;
	}

	residual.clear();
	p.clear(); 
	Ap.clear();
	apTT.clear();

	std::vector<double> cpuX(nSize);
	viennacl::copy(x, cpuX);
	x.clear();

	l->forEachLatticeSite([&](const int &i, const int &j, const int &k) {
		return cpuX[l->lIndex(i,j,k)];
	});
	return true;
}


template <class L, class GF, class I, class O>
bool fitPoints(
		const std::vector<vector3<I>> &samplePts, 
		const std::vector<O> &samples,
		const std::vector<std::tuple<int,int,int,O>> &optimizationVector,
		L *l,
		vector3<I> latticeShift,
		O *err = NULL, int *iter = NULL) {

	int nSize = l->numberOfLatticeSites();
	viennacl::compressed_matrix<O> phi(samplePts.size(), nSize), phiT(nSize, samplePts.size());
	createBigPhi<L,GF,I,O>(samplePts, l, latticeShift, phi, phiT);

	if(l->getLatticeName() == "bcc")
		ConjugateGradientBCC<L,GF,I,O>(phi, phiT, samples, optimizationVector, l, err, iter);
	else if(l->getLatticeName() == "cartesian")
		ConjugateGradientCC<L,GF,I,O>(phi, phiT, samples, optimizationVector, l, err, iter);
	else
		throw "fitPoints() - No cl kernel exists for this lattice type?";
	return true;
}

const char * _bcc_optimization_kernel_dbl =
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n"
"#define MODP(x, d) (((x)%(d) + (d)) % (d)) \n"
"__kernel void blnormbcc(\n"
"		  __global const double * vectorIn,\n"
"		  __global const int * xOffset, \n"
"		  __global const int * yOffset, \n"
"		  __global const int * zOffset, \n"
"		  __global const double * weights,\n"
"		  __global double * result,\n"
"		  unsigned int size, \n"
"		  unsigned int res, \n"
"		  double lambda ) { \n"
"		int res2 = res * 2; \n"

"	for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) { \n"
"		double temp = 0;	\n"

"		int z = (i % res2); \n"
"		int y = ((i - z)/res2) % res; \n"
"		int x = (((i - z)/res2) - y) / res; \n"
"		x = (x*2 + (z%2));\n"
"		y = (y*2 + (z%2));\n"

"		for(unsigned int j = 0; j < 175; j++) { \n"
"			int xx = (xOffset[j] + x); \n"
"			int yy = (yOffset[j] + y); \n"
"			int zz = (zOffset[j] + z); \n"
"			xx = (xx - abs(zz % 2))/ 2; \n"
"			yy = (yy - abs(zz % 2))/2; \n"

"			int index = zz + res2*(yy+res*xx); \n"

"			index = (xx < 0 || yy < 0 || zz < 0 || xx >= res2 || yy >= res2 || zz >= res2 || index >= size) ? -1 : index; \n"
"			temp += (index>=0)?vectorIn[index]*weights[j]:0; \n"
"		}\n"
"		result[i] += temp; \n"
"	}"
"};\n";

const char * _cc_optimization_kernel_dbl = 
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n"
"__kernel void blnorm(\n"
"       __global const double * vectorIn,\n"
"       __global const int * xOffset, \n"
"       __global const int * yOffset, \n"
"       __global const int * zOffset, \n"
"       __global const double * weights,\n"
"       __global double * result,\n"
"       unsigned int size, \n"
"		unsigned int nx, \n"
"		unsigned int ny, \n"
"		unsigned int nz, \n"
"		double lambda) { \n"
"	for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) { \n"
"		double temp = 0;	\n"
"		int z = i % nz; \n"
"		int y = ((i - z)/nz) % ny; \n"
"		int x = (((i - z)/nz) - y) / ny; \n"
"		for(unsigned int j = 0; j < 343; j++) { \n"
"			int xx = xOffset[j] + x; \n"
"			int yy = yOffset[j] + y; \n"
"			int zz = zOffset[j] + z; \n"
"			int index = zz + nx*(yy+ny*xx); \n"
"			if(xx < 0 || yy < 0 || zz < 0 || xx >= nx || yy >= ny || zz >= nz || index >= size) continue; \n"
"			temp += vectorIn[index]*weights[j]; \n"
"		}\n"
"		result[i] += temp; \n"
"	}"
"};\n";
};
};