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

// utility functions
#include <poisson/pointset.hpp>
#include <sisl/utility/isosurface.hpp>
#include <sisl/utility/scattered.hpp>

#include <tuple>

#define VESION_STRING "0.1"

using namespace sisl;
using namespace std;
using namespace TCLAP;

template <class T> void cartesianRecon(std::string input, std::string output, int res, bool shift, bool verbose, T lambda1, T lambda2, T scale, T stepsize);
template <class T> void bccRecon(std::string input, std::string output, int res, bool shift, bool verbose, T lambda1, T lambda2, T scale, T stepsize,  bool compat = true);
template <class T> void bcc4Recon(std::string input, std::string output, int res, bool shift, bool verbose, T lambda1, T lambda2, T scale, T stepsize, bool compat = true);
 
int main(int argc, char *argv[])
{
	try {
		CmdLine cmd("Poisson surface reconstruction", ' ', VESION_STRING);

		ValueArg<std::string> inputArg("i","input", "Input point file",true,"in.npts","filename");
		ValueArg<std::string> outputArg("o","output", "Output mesh/implicit file",true,"out.ply","filename");
		ValueArg<std::string> type("t","method", "Reconstruction type (CC, BCC, BCC4)",false,"BCC","string");
		ValueArg<int> res("r", "res", "Resolution of the reconstruction space", false, 256, "integer");
		ValueArg<double> lambda1("", "lambda", "Controls the compactness of the generating functions", false, 1.0, "float");
		ValueArg<double> lambda2("", "lambda2", "Controls the smoothness of the generating functions", false, 1.0, "float");
		ValueArg<double> scale("", "scale", "", false, 0.49, "float");
		ValueArg<double> ss("", "mstep", "", false, 0.005, "float");

		SwitchArg verboseSwitch("v","verbose","verbose console output.", cmd, false);
		SwitchArg shiftSwitch("s","shift","Shift reconstruction space.", cmd, false);
		SwitchArg compatSwitch("c","compat","fixes weighting to be compatible with an older version of the software that contained a bug.", cmd, false);

		cmd.add(ss);
		cmd.add(scale);
		cmd.add(lambda2);
		cmd.add(lambda1);
		cmd.add(res);
		cmd.add(type);
		cmd.add(outputArg);
		cmd.add(inputArg);

		cmd.parse(argc, argv);

		std::string in = inputArg.getValue();
		std::string out = outputArg.getValue();
		std::string method = type.getValue();
		int resolution = res.getValue();
		double l1 = lambda1.getValue();
		double l2 = lambda2.getValue();
		double scl = scale.getValue();
		double stepsize = ss.getValue();

		bool verbose = verboseSwitch.getValue();
		bool shift = shiftSwitch.getValue();
		bool compat = compatSwitch.getValue();

		if(method == "bcc")
			bccRecon<double>(in, out, resolution, shift, verbose, l1, l2, scl, stepsize, compat);
		else if (method == "bcc4")
			bcc4Recon<double>(in, out, resolution, shift, verbose, l1, l2, scl, stepsize, compat);
		else 
			cartesianRecon<double>(in, out, resolution, shift, verbose, l1, l2, scl, stepsize);

	}catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
	}catch (char const* e) {
		cerr << e << endl; 
	}
}

template <class T>
void dualConour(bcc_odd<T,T>){}

template <class T>
void cartesianRecon(std::string input, std::string output, int res, bool shift, bool verbose, T lambda1, T lambda2, T scale, T stepsize){
	typedef tp3cubic<T,T> GFType;
	typedef cartesian_odd<tp3cubic<T,T>, T, T> LATType;

	LATType *x, *y, *z, *l;
	utility::marchingCubes<T> surf;
	std::vector<std::tuple<int,int,int,T>> optimW, bl2, ac;

	T h = 1./T(res - 1);

	if(verbose) {
		cout << "Reconstruction on cartesian lattice. " << endl << "Input: " << input << endl;
		cout << "Output: " << output << endl << "Shift: " << shift << endl;
		cout << "Resoultion: " << res << endl;
		cout << "------" << endl << "Starting..." << endl;
	}

	// Read points
	if(verbose) cout << "Loading input pointset..." << endl;
	pointset<T> *p = new pointset<T>();
	if(!p->readXyz(input)) {
		throw "Couldn't read input file!";
	}

	vector3<T> c = p->getCenter();
	T radius = p->calculateRadius(c);
	matrix4x4<T> TS = matrix4x4<T>::translate(0.5, 0.5, 0.5) * matrix4x4<T>::uniformScale(scale/radius) * matrix4x4<T>::translate(-c);
	matrix4x4<T> invTS = matrix4x4<T>::translate(c) * matrix4x4<T>::uniformScale(radius/scale) * matrix4x4<T>::translate(-0.5, -0.5, -0.5);
	p->applyTransform(TS);
	printf("Loaded a point set with %d points.\n", p->count());

	// Fit points
	if(verbose) cout << "Optimizing gradient field... " << endl;
	x = new LATType(h);
	y = new LATType(h);
	z = new LATType(h);

	if(verbose) cout << "Creating mixed BL2-AC vector... " << endl;
	bl2 = GFType::getBeppoLevi2Norm();
	ac = GFType::autoCorrelation();
	for(unsigned int i = 0; i < ac.size(); i++){
		auto b = bl2[i];
		auto a = ac[i];

		if(std::get<0>(a) != std::get<0>(b) ||
			std::get<1>(a) != std::get<1>(b) ||
			std::get<2>(a) != std::get<2>(b))
			throw "Index mis-match while creating the optimization weight vector!";

		optimW.push_back(
			std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a),
				(std::get<3>(a)*lambda1*h + std::get<3>(b)*lambda2/(h*h))));
	}

	if(verbose) cout << "Creating big phi... " << endl;
	if(shift){
		int nSize = x->numberOfLatticeSites();
		viennacl::compressed_matrix<T> phi(p->getPositionVector().size(), nSize), phiT(nSize, p->getPositionVector().size());

		if(verbose) cout << "Solving systems... " << endl;

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0.5, 0, 0), phi, phiT);
		utility::ConjugateGradientCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(1, 0, 0)), 
			optimW, 
			x);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0, 0.5, 0), phi, phiT);
		utility::ConjugateGradientCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0, 1, 0)), 
			optimW, 
			y);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0, 0, 0.5), phi, phiT);
		utility::ConjugateGradientCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0, 0, 1)), 
			optimW, 
			z);

	}else{
		int nSize = x->numberOfLatticeSites();
		viennacl::compressed_matrix<T> phi(p->getPositionVector().size(), nSize), phiT(nSize, p->getPositionVector().size());
		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0, 0, 0), phi, phiT);

		if(verbose) cout << "Solving systems... " << endl;
		utility::ConjugateGradientCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(1, 0, 0)), 
			optimW, 
			x);

		utility::ConjugateGradientCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0, 1, 0)), 
			optimW, 
			y);

		utility::ConjugateGradientCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0, 0, 1)), 
			optimW, 
			z);
	}

	// Take Divergence
	if(verbose) cout << "Approximating divergence... " << endl;
	l = new LATType(h);

	if(shift)
		l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double sx = 1./(x->getScale()*12.);

			double dx = (
				 -1.*x->GV(i-1, j, k) + 
				  27.*x->GV(i, j, k) +
				 -27.*x->GV(i+1, j, k) + 
					 x->GV(i+2, j, k));
					
			double dy = (
				 (-1.*y->GV(i, j-1, k)) + 
				 ( 27.*y->GV(i, j, k)) + 
				 (-27.*y->GV(i, j+1, k)) + 
				 (	y->GV(i, j+2, k)));

			double dz = (
				 -1.*z->GV(i, j, k-1) + 
				  27.*z->GV(i, j, k) + 
				 -27.*z->GV(i, j, k+1) + 
				 z->GV(i, j, k+2)); 

			return (dx+dy+dz)*sx;
		});
	else
		l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double sx = 1./(x->getScale()*12.);

			double dx = (
			 -1.*x->GV(i-2, j, k) + 
			  8.*x->GV(i-1, j, k) +
			 -8.*x->GV(i+1, j, k) + 
				 x->GV(i+2, j, k)) * sx;
		
			double dy = (
			 (-1.*y->GV(i, j-2, k)) + 
			 ( 8.*y->GV(i, j-1, k)) + 
			 (-8.*y->GV(i, j+1, k)) + 
			 (	y->GV(i, j+2, k))) * sx;

			double dz = (
			 -1.*z->GV(i, j, k-2) + 
			  8.*z->GV(i, j, k-1) + 
			 -8.*z->GV(i, j, k+1) + 
				 z->GV(i, j, k+2)) * sx; 

			return dx+dy+dz;
		});

	// Cleanup
	delete x,y,z;

	// Invert Poisson
	if(verbose) cout << "Solving poisson... " << endl;
	l->frequencyFilter(tp3cubic<T,T>::poissonFilter);

	if(verbose) cout << "Guessing isovalue... " << endl;
	T isovalue = 0;
	std::vector<vector3<T>> s = p->getPositionVector();
	for(int i = 0; i < s.size(); i++)
		isovalue += l->f(s[i]);
	isovalue = isovalue / T(s.size());
	
	if(verbose) cout << "isovalue: " << isovalue << endl;

	delete p;

	if(verbose) cout << "Extracting isosurface... " << endl;
	
	surf.template marchLattice<LATType, T, T>(
		l, 
		NULL, 
		NULL, 
		NULL, 
		isovalue, 
		stepsize, 
		vector3<T>(0,0,0), 
		vector3<T>(1,1,1));

	// Transform to the input points
	surf.applyTransform(invTS);
	surf.writeSurface(output);
}

template <class T>
void bccRecon(std::string input, std::string output, int res, bool shift, bool verbose, T lambda1, T lambda2, T scale, T stepsize, bool compat = true){
	typedef quintic_box<T,T> GFType;
	typedef bcc_odd<quintic_box<T,T>, T, T> LATType;
	utility::marchingCubes<T> surf;
	LATType *x, *y, *z, *l;
	std::vector<std::tuple<int,int,int,T>> optimW, bl2, ac;

	T h = 1./(2*res);
	T compatH = 1./(2*res);

	if(verbose) {
		cout << "Reconstruction on BCC lattice. " << endl << "Input: " << input << endl;
		cout << "Output: " << output << endl << "Shift: " << shift << endl;
		cout << "Resoultion: " << res << endl;
	}

	// Read points
	if(verbose) cout << "Loading input pointset..." << endl;
	pointset<T> *p = new pointset<T>();
	if(!p->readXyz(input)) {
		throw "Couldn't read input file!";
	}

	vector3<T> c = p->getCenter();
	T radius = p->calculateRadius(c);
	matrix4x4<T> TS = matrix4x4<T>::translate(0.5, 0.5, 0.5) * matrix4x4<T>::uniformScale(scale/radius) * matrix4x4<T>::translate(-c);
	matrix4x4<T> invTS = matrix4x4<T>::translate(c) * matrix4x4<T>::uniformScale(radius/scale) * matrix4x4<T>::translate(-0.5, -0.5, -0.5);
	p->applyTransform(TS);
	printf("Loaded a point set with %d points.\n", p->count());

	// Fit points
	if(verbose) cout << "Optimizing gradient field... " << endl;
	x = new LATType(h);
	y = new LATType(h);
	z = new LATType(h);

	if(verbose) cout << "Creating mixed BL2-AC vector... " << endl;
	bl2 = quintic_box<T,T>::getBeppoLevi2Norm();
	ac = quintic_box<T,T>::autoCorrelation();
	for(unsigned int i = 0; i < ac.size(); i++){
		auto b = bl2[i];
		auto a = ac[i];
		if(std::get<0>(a) != std::get<0>(b) ||
			std::get<1>(a) != std::get<1>(b) ||
			std::get<2>(a) != std::get<2>(b))
			throw "Index mis-match while creating the optimization weight vector!";

		if(compat)
			optimW.push_back(
				std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a),
					(std::get<3>(a)*lambda1*compatH + std::get<3>(b)*lambda2/(compatH*compatH))
				)
			);
		else
			optimW.push_back(
				std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a),
					(std::get<3>(a)*lambda1*h + std::get<3>(b)*lambda2/(h*h))
				)
			);
	}

	if(verbose) cout << "Creating big phi... " << endl;
	if(shift){
		int nSize = x->numberOfLatticeSites();
		viennacl::compressed_matrix<T> phi(p->getPositionVector().size(), nSize), phiT(nSize, p->getPositionVector().size());
		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(-0.5,0.5,0.5), phi, phiT);
		
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0, 0.5, 0.5)), 
			optimW, 
			x);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0.5,-0.5,0.5), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.5, 0, 0.5)), 
			optimW, 
			y);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0.5,0.5,-0.5), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.5, 0.5, 0)), 
			optimW, 
			z);
	}else{
		int nSize = x->numberOfLatticeSites();
		viennacl::compressed_matrix<T> phi(p->getPositionVector().size(), nSize), phiT(nSize, p->getPositionVector().size());
		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0, 0, 0), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0, 0.5, 0.5)), 
			optimW, 
			x);

		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.5, 0, 0.5)), 
			optimW, 
			y);

		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.5, 0.5, 0)), 
			optimW, 
			z);
	}


	if(verbose) cout << "Taking divergence... " << endl;	

	l = new LATType(h);
	if(shift)
		l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double iscl = -1./(x->getScale()*12.); 

			double dx = 
				 -1.*x->GV(i-2, j+2, k+2) + 
				  27.*x->GV(i-1, j+1, k+1) +
				 -27.*x->GV(i, j, k) + 
					 x->GV(i+1, j-1, k-1);

			double dy = 
				 -1.*y->GV(i+2, j-2, k+2) + 
				  27.*y->GV(i+1, j-1, k+1) +
				 -27.*y->GV(i, j, k) + 
					 y->GV(i-1, j+1, k-1);
					
			double dz = 
				 -1.*z->GV(i+2, j+2, k-2) + 
				  27.*z->GV(i+1, j+1, k-1) +
				 -27.*z->GV(i, j, k) + 
					 z->GV(i-1, j-1, k+1);

			return iscl*(dx+dy+dz);
		});

	else
		l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double iscl = -1./(x->getScale()*12.); 

			double dx = 
				 -1.*x->GV(i-2, j+2, k+2) + 
				  8.*x->GV(i-1, j+1, k+1) +
				 -8.*x->GV(i+1, j-1, k-1) + 
					 x->GV(i+2, j-2, k-2);
			double dy = 
				 -1.*y->GV(i+2, j-2, k+2) + 
				  8.*y->GV(i+1, j-1, k+1) +
				 -8.*y->GV(i-1, j+1, k-1) + 
					 y->GV(i-2, j+2, k-2);
					
			double dz = 
				 -1.*z->GV(i+2, j+2, k-2) + 
				  8.*z->GV(i+1, j+1, k-1) +
				 -8.*z->GV(i-1, j-1, k+1) + 
					 z->GV(i-2, j-2, k+2);

			return iscl*(dx+dy+dz);
		});

	// Cleanup
	delete x,y,z;

	// Invert Poisson
	if(verbose) cout << "Solving poisson... " << endl;
	l->frequencyFilter(GFType::poissonFilter);

	if(verbose) cout << "Guessing isovalue ... " << endl;
	T isovalue = 0;
	std::vector<vector3<T>> s = p->getPositionVector();
	for(int i = 0; i < s.size(); i++)
		isovalue += l->f(s[i]);
	isovalue = isovalue / T(s.size());
	
	if(verbose) cout << "isovalue: " << isovalue << endl;

	delete p;

	if(verbose) cout << "Extracting isosurface... " << endl;
	
	surf.template marchLattice<LATType, T, T>(
		l, 
		NULL, 
		NULL, 
		NULL, 
		isovalue, 
		stepsize, 
		vector3<T>(0,0,0), 
		vector3<T>(1,1,1));

	// Transform to the input points
	surf.applyTransform(invTS);
	surf.writeSurface(output);
}

template <class T>
void bcc4Recon(std::string input, std::string output, int res, bool shift, bool verbose, T lambda1, T lambda2, T scale, T stepsize, bool compat = true){
	typedef quintic_box<T,T> GFType;
	typedef bcc_odd<quintic_box<T,T>, T, T> LATType;
	utility::marchingCubes<T> surf;
	LATType *x, *y, *z, *w, *l;
	std::vector<std::tuple<int,int,int,T>> optimW, bl2, ac;

	T h = 1./(2*res);
	T compatH = 1./(2*res);

	if(verbose) {
		cout << "Reconstruction on BCC lattice. " << endl << "Input: " << input << endl;
		cout << "Output: " << output << endl << "Shift: " << shift << endl;
		cout << "Resoultion: " << res << endl;
	}

	// Read points
	if(verbose) cout << "Loading input pointset..." << endl;
	pointset<T> *p = new pointset<T>();
	if(!p->readXyz(input)) {
		throw "Couldn't read input file!";
	}

	vector3<T> c = p->getCenter();
	T radius = p->calculateRadius(c);
	matrix4x4<T> TS = matrix4x4<T>::translate(0.5, 0.5, 0.5) * matrix4x4<T>::uniformScale(scale/radius) * matrix4x4<T>::translate(-c);
	matrix4x4<T> invTS = matrix4x4<T>::translate(c) * matrix4x4<T>::uniformScale(radius/scale) * matrix4x4<T>::translate(-0.5, -0.5, -0.5);
	p->applyTransform(TS);
	printf("Loaded a point set with %d points.\n", p->count());

	// Fit points
	if(verbose) cout << "Optimizing gradient field... " << endl;
	x = new LATType(h);
	y = new LATType(h);
	z = new LATType(h);
	w = new LATType(h);

	if(verbose) cout << "Creating mixed BL2-AC vector... " << endl;
	bl2 = quintic_box<T,T>::getBeppoLevi2Norm();
	ac = quintic_box<T,T>::autoCorrelation();
	for(unsigned int i = 0; i < ac.size(); i++){
		auto b = bl2[i];
		auto a = ac[i];

		if(std::get<0>(a) != std::get<0>(b) ||
			std::get<1>(a) != std::get<1>(b) ||
			std::get<2>(a) != std::get<2>(b))
			throw "Index mis-match while creating the optimization weight vector!";

		if(compat)
			optimW.push_back(
				std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a),
					(std::get<3>(a)*lambda1*compatH + std::get<3>(b)*lambda2/(compatH*compatH))
				)
			);
		else
			optimW.push_back(
				std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a),
					(std::get<3>(a)*lambda1*h + std::get<3>(b)*lambda2/(h*h))
				)
			);
	}

	if(verbose) cout << "Creating big phi... " << endl;
	if(shift){
		int nSize = x->numberOfLatticeSites();
		viennacl::compressed_matrix<T> phi(p->getPositionVector().size(), nSize), phiT(nSize, p->getPositionVector().size());
		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(-0.5,0.5,0.5), phi, phiT);
		
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(-0.25, 0.25, 0.25)), 
			optimW, 
			x);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0.5,-0.5,0.5), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.25, -0.25, 0.25)), 
			optimW, 
			y);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0.5,0.5,-0.5), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.25, 0.25, -0.25)), 
			optimW, 
			z);

		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(-0.5,-0.5,-0.5), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(-0.25, -0.25, -0.25)), 
			optimW, 
			w);
	}else{
		int nSize = x->numberOfLatticeSites();
		viennacl::compressed_matrix<T> phi(p->getPositionVector().size(), nSize), phiT(nSize, p->getPositionVector().size());
		utility::createBigPhi<LATType,GFType,T,T>(p->getPositionVector(), x, sisl::vector3<T>(0, 0, 0), phi, phiT);
		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(-0.25, 0.25, 0.25)), 
			optimW, 
			x);

		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT,
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.25, -0.25, 0.25)), 
			optimW, 
			y);

		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(0.25, 0.25, -0.25)), 
			optimW, 
			z);

		utility::ConjugateGradientBCC<LATType,GFType,T,T>(
			phi, phiT, 
			p->getWeigthedNormalColumnVector(sisl::vector3<T>(-0.25, -0.25, -0.25)), 
			optimW, 
			w);
	}


	if(verbose) cout << "Taking divergence... " << endl;	

	l = new LATType(h);
	if(shift)
		l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double iscl = -1./(x->getScale()*12.); 

			double dx = 
				 -1.*x->GV(i-2, j+2, k+2) + 
				  27.*x->GV(i-1, j+1, k+1) +
				 -27.*x->GV(i, j, k) + 
					 x->GV(i+1, j-1, k-1);

			double dy = 
				 -1.*y->GV(i+2, j-2, k+2) + 
				  27.*y->GV(i+1, j-1, k+1) +
				 -27.*y->GV(i, j, k) + 
					 y->GV(i-1, j+1, k-1);
					
			double dz = 
				 -1.*z->GV(i+2, j+2, k-2) + 
				  27.*z->GV(i+1, j+1, k-1) +
				 -27.*z->GV(i, j, k) + 
					 z->GV(i-1, j-1, k+1);


			double dw = 
				 -1.*w->GV(i-2, j-2, k-2) + 
				  27.*w->GV(i-1, j-1, k-1) +
				 -27.*w->GV(i, j, k) + 
					 w->GV(i+1, j+1, k+1);

			return iscl*(dx+dy+dz+dw);
		});

	else
		l->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
			double iscl = -1./(x->getScale()*12.); 

			double dx = 
				 -1.*x->GV(i-2, j+2, k+2) + 
				  8.*x->GV(i-1, j+1, k+1) +
				 -8.*x->GV(i+1, j-1, k-1) + 
					 x->GV(i+2, j-2, k-2);
			double dy = 
				 -1.*y->GV(i+2, j-2, k+2) + 
				  8.*y->GV(i+1, j-1, k+1) +
				 -8.*y->GV(i-1, j+1, k-1) + 
					 y->GV(i-2, j+2, k-2);
					
			double dz = 
				 -1.*z->GV(i+2, j+2, k-2) + 
				  8.*z->GV(i+1, j+1, k-1) +
				 -8.*z->GV(i-1, j-1, k+1) + 
					 z->GV(i-2, j-2, k+2);

			double dw = 
				 -1.*w->GV(i-2, j-2, k-2) + 
				  8.*w->GV(i-1, j-1, k-1) +
				 -8.*w->GV(i+1, j+1, k+1) + 
					 w->GV(i+2, j+2, k+2);
			return iscl*(dx+dy+dz+dw);
		});

	// Cleanup
	delete x,y,z,w;

	// Invert Poisson
	if(verbose) cout << "Solving poisson... " << endl;
	l->frequencyFilter(GFType::poissonFilter);

	if(verbose) cout << "Guessing isovalue ... " << endl;
	T isovalue = 0;
	std::vector<vector3<T>> s = p->getPositionVector();
	for(int i = 0; i < s.size(); i++)
		isovalue += l->f(s[i]);
	isovalue = isovalue / T(s.size());
	
	if(verbose) cout << "isovalue: " << isovalue << endl;

	delete p;

	if(verbose) cout << "Extracting isosurface... " << endl;
	
	surf.template marchLattice<LATType, T, T>(
		l, 
		NULL, 
		NULL, 
		NULL, 
		isovalue, 
		stepsize, 
		vector3<T>(0,0,0), 
		vector3<T>(1,1,1));

	// Transform to the input points
	surf.applyTransform(invTS);
	surf.writeSurface(output);
}
