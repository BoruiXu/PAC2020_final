#include "common.h"

TDEF(fftw)
TDEF(nufft)

/* Constructor */
NUFFT3D::NUFFT3D(int N,int OF,float* wx,float* wy,float* wz,int P,int prechopX,int prechopY,int prechopZ,int postchopX,int postchopY,int postchopZ,int offsetX,int offsetY,int offsetZ,int W,int L) {

	// Assignments
	this->N = N;
	this->OF = OF;
	N2 = N*OF;
	this->wx = wx;
	this->wy = wy;
	this->wz = wz;
	this->P = P;
	this->prechopX = prechopX;
	this->prechopY = prechopY;
	this->prechopZ = prechopZ;
	this->postchopX = postchopX;
	this->postchopY = postchopY;
	this->postchopZ = postchopZ;
	this->offsetX = offsetX;
	this->offsetY = offsetY;
	this->offsetZ = offsetZ;
	this->W = W;
	this->L = L;


	my_desc1_handle = NULL;
	my_desc2_handle = NULL;
	
	f = (complex<float>*)memalign(16,N2*N2*N2*sizeof(complex<float>));

	
	MKL_LONG DIMS[3] = {N2,N2,N2};

    status = DftiCreateDescriptor(&my_desc1_handle, DFTI_SINGLE,
                                  DFTI_COMPLEX, 3, DIMS);
    status = DftiCommitDescriptor(my_desc1_handle);

    status = DftiCreateDescriptor(&my_desc2_handle, DFTI_SINGLE,
                                  DFTI_COMPLEX, 3, DIMS);
    status = DftiCommitDescriptor(my_desc2_handle);

	//fwdPlan = fftwf_plan_dft(3,DIMS,reinterpret_cast<fftwf_complex*>(f),reinterpret_cast<fftwf_complex*>(f),FFTW_FORWARD,FFTW_ESTIMATE);
	//adjPlan = fftwf_plan_dft(3,DIMS,reinterpret_cast<fftwf_complex*>(f),reinterpret_cast<fftwf_complex*>(f),FFTW_BACKWARD,FFTW_ESTIMATE);
	buildLUT();
	getScalingFunction();
}

/* Destructor */
NUFFT3D::~NUFFT3D() {
	//fftwf_destroy_plan(fwdPlan);
	//fftwf_destroy_plan(adjPlan);
	status = DftiFreeDescriptor(&my_desc1_handle);
	status = DftiFreeDescriptor(&my_desc2_handle);
	free(f);
	f = NULL;
}

/* Initialize multithreaded FFTW (p threads) */
void NUFFT3D::init(int nThreads) {
	//fftwf_init_threads();
	//fftwf_plan_with_nthreads(nThreads);
}

/* Forward NUFFT transform */
void NUFFT3D::fwd(complex<float>* u,complex<float>* raw) {
	// Apodization and zero-padding
	int startX = N*(OF-1)/2;
	int startY = N*(OF-1)/2;
	int startZ = N*(OF-1)/2;
	TSTART(nufft);
	TSTART(fftw);
	for (int i=0;i<N2*N2*N2;i++) {f[i] = 0;}
	TEND(fftw); TPRINT(fftw, "  Init_F FWD");
	TSTART(fftw);
	//加权  q在getScalingFunction();中获得
	for (int x=0;x<N;x++) {
		for (int y=0;y<N;y++) {
			for (int z=0;z<N;z++) {
				f[(x+startX+offsetX)*N2*N2+(y+startY+offsetY)*N2+(z+startZ+offsetZ)] = u[x*N*N+y*N+z]/q[x*N*N+y*N+z];
			}
		}
	}
	TEND(fftw); TPRINT(fftw, "  Roundoff_Corr FWD");

	// (Oversampled) FFT    过采样: f
	TSTART(fftw);
	chop3D(f,N2,N2,N2,prechopX,prechopY,prechopZ);
	TEND(fftw); TPRINT(fftw, "  PreChop FWD")
	TSTART(fftw);
	//fftwf_execute(fwdPlan);
	status = DftiComputeForward(my_desc1_handle, f);
	TEND(fftw); TPRINT(fftw, "  FFTW FWD");
	TSTART(fftw)
	chop3D(f,N2,N2,N2,postchopX,postchopY,postchopZ);
	TEND(fftw); TPRINT(fftw, "  PostChop FWD");

	// Pull from grid
	TSTART(fftw);
	for (int p=0;p<P;p++) {

		int kx2[2*W+1];
		int ky2[2*W+1];
		int kz2[2*W+1];
		float winX[2*W+1];
		float winY[2*W+1];
		float winZ[2*W+1];

		// Form x interpolation kernel
		float kx = N2*(wx[p]+0.5); 
		int x1 = (int)ceil(kx-W); 
		int x2 = (int)floor(kx+W); 
		int lx = x2-x1+1;
		for (int nx=0;nx<lx;nx++) {
			kx2[nx] = mod(nx+x1,N2);
			winX[nx] = LUT[(int)round(((L-1)/W)*abs(nx+x1-kx))];
		}

		// Form y interpolation kernel
		float ky = N2*(wy[p]+0.5);
		int y1 = (int)ceil(ky-W);
		int y2 = (int)floor(ky+W);
		int ly = y2-y1+1;
		for (int ny=0;ny<ly;ny++) {
			ky2[ny] = mod(ny+y1,N2);
			winY[ny] = LUT[(int)round(((L-1)/W)*abs(ny+y1-ky))];
		}

		// Form z interpolation kernel
		float kz = N2*(wz[p]+0.5);
		int z1 = (int)ceil(kz-W);
		int z2 = (int)floor(kz+W);
		int lz = z2-z1+1;
		for (int nz=0;nz<lz;nz++) {
			kz2[nz] = mod(nz+z1,N2);
			winZ[nz] = LUT[(int)round(((L-1)/W)*abs(nz+z1-kz))];
		}

		// Interpolation
		raw[p] = 0;
		for (int nx=0;nx<lx;nx++) {
			for (int ny=0;ny<ly;ny++) {
				for (int nz=0;nz<lz;nz++) {
					raw[p] += f[kx2[nx]*N2*N2+ky2[ny]*N2+kz2[nz]]*winX[nx]*winY[ny]*winZ[nz];
				}
			}
		}
	}
	TEND(fftw); TPRINT(fftw, "  Convolution FWD");
	TEND(nufft); TPRINT(nufft, "NUFFT FWD");
}

/* Adjoint NUFFT transform */
void NUFFT3D::adj(complex<float>* raw,complex<float>* u) {
	TSTART(nufft);

	// Push to grid
	TSTART(fftw);
	for (int i=0;i<N2*N2*N2;i++) {f[i] = 0;}
	TEND(fftw); TPRINT(fftw, "  Init_F ADJ");
	TSTART(fftw)
	for (int p=0;p<P;p++) {

		int kx2[2*W+1];
		int ky2[2*W+1];
		int kz2[2*W+1];
		float winX[2*W+1];
		float winY[2*W+1];
		float winZ[2*W+1];

		// Form x interpolation kernel
		float kx = N2*(wx[p]+0.5);
		int x1 = (int)ceil(kx-W);
		int x2 = (int)floor(kx+W);
		int lx = x2-x1+1;
		for (int nx=0;nx<lx;nx++) {
			kx2[nx] = mod(nx+x1,N2);
			winX[nx] = LUT[(int)round(((L-1)/W)*abs(nx+x1-kx))];
		}

		// Form y interpolation kernel
		float ky = N2*(wy[p]+0.5);
		int y1 = (int)ceil(ky-W);
		int y2 = (int)floor(ky+W);
		int ly = y2-y1+1;
		for (int ny=0;ny<ly;ny++) {
			ky2[ny] = mod(ny+y1,N2);
			winY[ny] = LUT[(int)round(((L-1)/W)*abs(ny+y1-ky))];
		}

		// Form z interpolation kernel
		float kz = N2*(wz[p]+0.5);
		int z1 = (int)ceil(kz-W);
		int z2 = (int)floor(kz+W);
		int lz = z2-z1+1;
		for (int nz=0;nz<lz;nz++) {
			kz2[nz] = mod(nz+z1,N2);
			winZ[nz] = LUT[(int)round(((L-1)/W)*abs(nz+z1-kz))];
		}

		// Interpolation
		for (int nx=0;nx<lx;nx++) {
			for (int ny=0;ny<ly;ny++) {
				for (int nz=0;nz<lz;nz++) {
					f[kx2[nx]*N2*N2+ky2[ny]*N2+kz2[nz]] += raw[p]*winX[nx]*winY[ny]*winZ[nz];
				}
			}
		}
	}
	TEND(fftw); TPRINT(fftw, "  Convolution ADJ");
	// (Oversampled) FFT
	TSTART(fftw);
	chop3D(f,N2,N2,N2,postchopX,postchopY,postchopZ);
	TEND(fftw); TPRINT(fftw, "  PostChop ADJ");
	TSTART(fftw);
	//fftwf_execute(adjPlan);
	status = DftiComputeBackward(my_desc2_handle, f);	
	TEND(fftw); TPRINT(fftw, "  FFTW ADJ");
	TSTART(fftw);
	chop3D(f,N2,N2,N2,prechopX,prechopY,prechopZ);
	TEND(fftw); TPRINT(fftw, "  PreChop ADJ");
	// Deapodize and truncate
	int startX = N*(OF-1)/2;
	int startY = N*(OF-1)/2;
	int startZ = N*(OF-1)/2;
	TSTART(fftw);
	for (int i=0;i<N*N*N;i++) {u[i] = 0;}
	TEND(fftw); TPRINT(fftw, "  Init_U ADJ");
	TSTART(fftw)
	for (int x=0;x<N;x++) {
		for (int y=0;y<N;y++) {
			for (int z=0;z<N;z++) {
				u[x*N*N+y*N+z] = f[(x+startX+offsetX)*N2*N2+(y+startY+offsetY)*N2+(z+startZ+offsetZ)]/q[x*N*N+y*N+z];
			}
		}
	}
	TEND(fftw); TPRINT(fftw, "  Roundoff_Corr ADJ");
	TEND(nufft); TPRINT(nufft, "NUFFT ADJ")
	return;
}

/* Internal lookup table generation function for interpolation kernel (Kaiser-Bessel) */
void NUFFT3D::buildLUT() {
	LUT = new float[L];
	float *d = linspace<float>(0,W,L);
	float pi = boost::math::constants::pi<float>();
	float alpha = pi*sqrt(((2*(float)W/OF)*(OF-0.5))*((2*(float)W/OF)*(OF-0.5))-0.8);
	for (int l=0;l<L;l++) {
		LUT[l] = boost::math::cyl_bessel_i(0,alpha*sqrt(1-(d[l]*d[l])/(W*W))) / boost::math::cyl_bessel_i(0,alpha);
	}
}

/* Internal scaling generation function */
void NUFFT3D::getScalingFunction() {

	float dx,dy,dz;
	float s = 0;

	// Create a volume with a copy of the interpolation kernel centered at the origin, then normalize
	for (int i=0;i<N2*N2*N2;i++) {f[i] = 0;}
	for (int x=N2/2-W;x<=N2/2+W;x++) {
        dx = abs(((float)x-N2/2)/W); 
		for (int y=N2/2-W;y<=N2/2+W;y++) {
			dy = abs(((float)y-N2/2)/W);
			for (int z=N2/2-W;z<=N2/2+W;z++) {
				dz = abs(((float)z-N2/2)/W);
				f[x*N2*N2+y*N2+z] = complex<float>(LUT[(int)round((L-1)*dx)]*LUT[(int)round((L-1)*dy)]*LUT[(int)round((L-1)*dz)],0);
				s = s + norm(f[x*N2*N2+y*N2+z]);
			}
		}
	}
	s = sqrt(s);
	for (int x=N2/2-W;x<=N2/2+W;x++) {
		for (int y=N2/2-W;y<=N2/2+W;y++) {
			for (int z=N2/2-W;z<=N2/2+W;z++) {
				f[x*N2*N2+y*N2+z] = f[x*N2*N2+y*N2+z]/s;
			}
		}
	}
	
	// (Oversampled) FFT
	chop3D(f,N2,N2,N2,postchopX,postchopY,postchopZ);
	//fftwf_execute(adjPlan);
	status = DftiComputeBackward(my_desc2_handle, f);	
	chop3D(f,N2,N2,N2,prechopX,prechopY,prechopZ);

	// Truncate and keep only the real component (presuming Fourier domain symmetry)
	q = new float[N*N*N];
	int startX = N*(OF-1)/2;
	int startY = N*(OF-1)/2;
	int startZ = N*(OF-1)/2;
	for (int i=0;i<N*N*N;i++) {q[i] = 0;}
	for (int x=0;x<N;x++) {
		for (int y=0;y<N;y++) {
			for (int z=0;z<N;z++) {
				q[x*N*N+y*N+z] = real(f[(x+startX+offsetX)*N2*N2+(y+startY+offsetY)*N2+(z+startZ+offsetZ)]);
			}
		}
	}

	return;
}
