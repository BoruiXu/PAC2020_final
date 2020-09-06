#include "common.h"

chrono::high_resolution_clock::time_point fwd_t0, fwd_t1, adj_t0, adj_t1, wall_t0, wall_t1;

using namespace std;

int main(int argc, char **argv) {
	char rawfn[256], coordfn[256], ksfn[256], reffn[256], outfn[256];
	char base[256];
	int K, S;
	int numIters = 20;

	int numThreads = 1;
	if(argc < 2) {
          printf("Usages: %s <basename> [numThreads] [numIterations]\n", argv[0]);
        }
	strncpy(base, argv[1], 255);
	snprintf(ksfn, 256, "%s.%s", base, "KS");
	snprintf(rawfn, 256, "%s.%s", base, "raw");
	snprintf(coordfn, 256, "%s.%s", base, "coord");

	if(argc > 2) {
		numThreads = atoi(argv[2]);
		if(numThreads > 0) {
			omp_set_num_threads(numThreads);
		}
	}
	numThreads = omp_get_max_threads();
	printf("Using %d Threads\n", numThreads);
	if(argc > 3) {
		numIters = atoi(argv[3]);
		if (numIters > 20) numIters = 20;
		if (numIters < 1) printf("Number of interarions must be > 0.\n");
	}
	/* ****************** */
	/* *** Load files *** */
	/* ********************************************************************* */

	fstream in;

	// Read raw data file
	in.open(ksfn, ios::in);
	if(!in) {
		printf("Can't open file \"%s\"\n", ksfn);
		exit(1);
	}
	in >> K >> S;
	in.close();
	printf("K = %d , S = %d, TotalSamples = %d\n", K, S, K*S);
	float *temp = new float[K*S*2];
	
	in.open(rawfn,ios::in|ios::binary);
	if(!in) {
		printf("Can't open file \"%s\"\n", rawfn);
		exit(1);
	}
	in.seekg(0,ios::beg);	
	in.read(reinterpret_cast<char*>(temp),K*S*2*sizeof(float));
	complex<float> *raw = new complex<float>[K*S];
	for (int s=0;s<S;s++) {
		for (int k=0;k<K;k++) {
			raw[s*K+k] = complex<float>(temp[2*(s*K+k)],temp[2*(s*K+k)+1]);
		}
	}
	delete[] temp; temp = NULL;
	in.close();

	// Read sampling file	
	float *temp2 = new float[K*S*3];
	in.open(coordfn,ios::in|ios::binary);
	if(!in) {
		printf("Can't open file \"%s\"\n", coordfn);
		exit(1);
	}
	in.seekg(0,ios::beg);
	in.read(reinterpret_cast<char*>(temp2),K*S*3*sizeof(float));
	float *wx = new float[K*S]; float *wy = new float[K*S];	float *wz = new float[K*S];
	for (int s=0;s<S;s++) {
		for (int k=0;k<K;k++) {
			wx[s*K+k] = (0.1)*temp2[0*S*K+s*K+k];
			wy[s*K+k] = (0.1)*temp2[1*S*K+s*K+k];
			wz[s*K+k] = (0.1)*temp2[2*S*K+s*K+k];
		}
	}
	delete[] temp2; temp2 = NULL;
	in.close();

	/* ********************************* */
	/* *** Initialize NUFFT Operator *** */
	/* ********************************************************************* */

	// Multithreaded initialization
	NUFFT3D::init(numThreads);

	// Initialize NUFFT operator specifics
	int N = 220;
	int W = 4;
	int OF = 2;
	int L = 25000;
	int P = K*S;
	int prechopX = 1;
	int postchopX = 1;
	int prechopY = 1;
	int postchopY = 1;
	int prechopZ = 1;
	int postchopZ = 1;
	int offsetX = 0;
	int offsetY = 0;
	int offsetZ = -35;
	NUFFT3D *nufft = new NUFFT3D(N,OF,wx,wy,wz,P,prechopX,prechopY,prechopZ,postchopX,postchopY,postchopZ,offsetX,offsetY,offsetZ,W,L);

	/* ******************************** */
	/* *** Begin CG Iterative Recon *** */
	/* ********************************************************************* */	

	float lambda = 1e-10;
	int maxCGIter = numIters;
	float epsilon = 1e-30;

	// Standard CG algorithm for solving (A'A+lambdaI)u = A'f;
	complex<float> alpha,den;
	float beta,delta_old,delta;
	complex<float> *u = new complex<float>[N*N*N];
	complex<float> *r = new complex<float>[N*N*N];
	complex<float> *p = new complex<float>[N*N*N];
	complex<float> *Ap = new complex<float>[N*N*N];
	complex<float> *z = new complex<float>[P];

	wall_t0=Clock::now();
	adj_t0=Clock::now();
	nufft->adj(raw,r);
	adj_t1=Clock::now();
	cout << "ADJ time = \t\t" << chrono::duration_cast<chrono::microseconds>(adj_t1 - adj_t0).count()/1e6 << "\t\t sec" << endl;
	cout << "=================================================" << endl;
	for (int i=0;i<N*N*N;i++) {
		u[i] = 0;
		p[i] = r[i];
	}
	delta_old = 0; for (int i=0;i<N*N*N;i++) delta_old += norm(r[i]);
	for (int iter=0;iter<maxCGIter;iter++) {
		cout << "Iteration " << iter+1 << " :" << endl;
		fwd_t0=Clock::now();
		nufft->fwd(p,z);
		fwd_t1=Clock::now();
		cout << "FWD time = \t\t" << chrono::duration_cast<chrono::microseconds>(fwd_t1 - fwd_t0).count()/1e6 << "\t\t sec" << endl;
		cout << "-------------------------------------------------" << endl;
		adj_t0=Clock::now();
		nufft->adj(z,Ap);
		adj_t1=Clock::now();
		cout << "ADJ time = \t\t" << chrono::duration_cast<chrono::microseconds>(adj_t1 - adj_t0).count()/1e6 << "\t\t sec" << endl;
		for (int i=0;i<N*N*N;i++) Ap[i] += lambda*p[i];
		den = epsilon; for (int i=0;i<N*N*N;i++) den += conj(p[i])*Ap[i];
		alpha = delta_old/den;
		for (int i=0;i<N*N*N;i++) {
			u[i] = u[i] + alpha*p[i];
			r[i] = r[i] - alpha*Ap[i];
		}
		delta = 0; for (int i=0;i<N*N*N;i++) delta += norm(r[i]);
		beta = delta/(delta_old+epsilon);
		delta_old = delta;
		for (int i=0;i<N*N*N;i++) p[i] = r[i] + beta*p[i];
		cout << "Iteration " << iter+1 << " out of " << maxCGIter << endl;
		cout << "=================================================" << endl;
	}
	wall_t1=Clock::now();
	cout << "Total Computation Time:" << endl << "\t\t" << chrono::duration_cast<chrono::microseconds>(wall_t1 - wall_t0).count()/1e6 << "\t sec" << endl;
	delete nufft;

	/* ******************** */
	/* *** Save results *** */
	/* ********************************************************************* */
	snprintf(outfn, 256, "%s.%s.%d", base, "bin", maxCGIter);
	fstream out;
	out.open(outfn, ios::out|ios::binary);
	out.seekp(0,ios::beg);
	out.write(reinterpret_cast<char*>(u),N*N*N*sizeof(complex<float>));
	out.close();

	return 0;
}
