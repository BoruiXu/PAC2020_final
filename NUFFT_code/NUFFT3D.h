#ifndef NUFFT3D_H
#define NUFFT3D_H

#include "MathOps.h"
#include <complex>
#include <fftw3.h>
#include <malloc.h>

using namespace std;

class NUFFT3D {
	private:
		complex<float> *f;
		int N;
		int OF;
		int N2;
		float* wx;
		float* wy;
		float* wz;
		int P;
		int prechopX;
		int prechopY;
		int prechopZ;
		int postchopX;
		int postchopY;
		int postchopZ;
		int offsetX;
		int offsetY;
		int offsetZ;
		int W;
		int L;
		fftwf_plan fwdPlan;
		fftwf_plan adjPlan;
		float* LUT;
		void buildLUT();
		void getScalingFunction();
		float *q;
	public:
		NUFFT3D(int,int,float*,float*,float*,int,int,int,int,int,int,int,int,int,int,int,int);
		~NUFFT3D();
		static void init(int);
		void fwd(complex<float>*,complex<float>*);
		void adj(complex<float>*,complex<float>*);
};

#endif