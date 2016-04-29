#ifndef LC_HPP
#define LC_HPP

#include <complex>
#include <mkl_types.h>
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

//#include "CARMA.hpp"
//#include "MCMC.hpp"
//#include "Constants.hpp"
//#include "rdrand.hpp"

using namespace std;

	class LCData {
	public:
		int numCadences;
		double dt;
		double tolIR;
		double fracIntrinsicVar;
		double fracNoiseToSignal;
		double maxSigma;
		double minTimescale;
		double maxTimescale;
		double *t;
		double *x;
		double *y;
		double *yerr;
		double *mask;
		double *lcX;
		double *lcP;
		LCData();
		~LCData();
	};

#endif