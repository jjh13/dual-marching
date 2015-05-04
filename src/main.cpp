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

#include <tuple>

#define VESION_STRING "0.1"

using namespace sisl;
using namespace std;
using namespace TCLAP;

template <class T>
void dualConour(bcc_odd<quintic_box<T,T>, T, T>  *lattice);

int main(int argc, char *argv[])
{
	try {
		CmdLine cmd("Dual BCC marching", ' ', VESION_STRING);

		ValueArg<std::string> outputArg("o", "output", "Output mesh/implicit file",true,"out.ply","filename");

		// ValueArg<std::string> type("t","method", "Reconstruction type (CC, BCC, BCC4)",false,"BCC","string");
		// ValueArg<int> res("r", "res", "Resolution of the reconstruction space", false, 256, "integer");
		// ValueArg<double> lambda1("", "lambda", "Controls the compactness of the generating functions", false, 1.0, "float");
		// ValueArg<double> lambda2("", "lambda2", "Controls the smoothness of the generating functions", false, 1.0, "float");
		// ValueArg<double> scale("", "scale", "", false, 0.49, "float");
		// ValueArg<double> ss("", "mstep", "", false, 0.005, "float");

		// SwitchArg verboseSwitch("v","verbose","verbose console output.", cmd, false);
		// SwitchArg shiftSwitch("s","shift","Shift reconstruction space.", cmd, false);
		// SwitchArg compatSwitch("c","compat","fixes weighting to be compatible with an older version of the software that contained a bug.", cmd, false);

		// cmd.add(ss);
		// cmd.add(scale);
		// cmd.add(lambda2);
		// cmd.add(lambda1);
		// cmd.add(res);
		// cmd.add(type);
		cmd.add(outputArg);
		// cmd.add(inputArg);

		cmd.parse(argc, argv);

		// std::string in = inputArg.getValue();
		std::string output = outputArg.getValue();
		// std::string method = type.getValue();
		// int resolution = res.getValue();
		// double l1 = lambda1.getValue();
		// double l2 = lambda2.getValue();
		// double scl = scale.getValue();
		// double stepsize = ss.getValue();

		// bool verbose = verboseSwitch.getValue();
		// bool shift = shiftSwitch.getValue();
		// bool compat = compatSwitch.getValue();

		cout << output << endl;

		typedef quintic_box<double, double> GFType;
		typedef bcc_odd<quintic_box<double, double>, double, double> LATType;
		LATType *lat = new LATType(1/128);

		dualConour(lat);


	}catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
	}catch (char const* e) {
		cerr << e << endl; 
	}
}

template <class T>
void dualConour(bcc_odd<quintic_box<T,T>, T, T>  *lattice){

	// typedef linear_bcc_box<T,T> GFType;
	// typedef bcc_odd<linear_bcc_box<T,T>, T, T> LATType;

	/**********************
	 * 
	 **********************/
	lattice->forEachLatticeSite([&](const int & i, const int &j, const int &k) { 
		printf("%d, %d, %d\n", i, j, k);
	});
	return;
}


