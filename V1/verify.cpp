#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <complex>

using namespace std;

int N = 220;

int main(int argc, char **argv) {
	complex<float> *u = new complex<float>[N*N*N];
	complex<float> *r = new complex<float>[N*N*N];
        int pass = 1;
        int errors = 0;
        complex<float> avg = 0, avgref = 0;

	if(argc < 2 || argc > 3) {
          printf("Usages: %s <reference> <result>\n", argv[0]);
          exit(1);
        }

	/* ****************** */
	/* *** Load files *** */
	/* ********************************************************************* */

	fstream refin, resultin;
	printf("reference file: %s\n", argv[1]);
        refin.open(argv[1], ios::in|ios::binary);
        if(!refin) {
                printf("Can't open reference file.\n");
                exit(2);
        }
        refin.seekp(0,ios::beg);
        refin.read(reinterpret_cast<char*>(r),N*N*N*sizeof(complex<float>));
        refin.close();

	printf("result file: $s\n", argv[2]);
        resultin.open(argv[2], ios::in|ios::binary);
        if(!resultin) {
                printf("Can't open result file.\n");
                exit(3);
        }
        resultin.seekp(0,ios::beg);
        resultin.read(reinterpret_cast<char*>(u),N*N*N*sizeof(complex<float>));
        resultin.close();

	/* ******************** */
	/* *** Verify results *** */
	/* ********************************************************************* */
        for(int i = 0; i < N*N*N; i++) {
		avg += u[i];
		avgref += r[i];
		if(isnan(real(u[i])) || isnan(imag(u[i])) || (norm(r[i] - u[i])/norm(r[i]) > 1.0e-3)) {
			errors++;
			if(errors <= 10) {
				printf("Error #%03d @ %d = %g\n", errors, i, norm(r[i] - u[i])/norm(r[i]));
				printf("ref[%d] = (%g,%g),  computed[%d] = (%g,%g)\n", i, real(r[i]), imag(r[i]), i, real(u[i]), imag(u[i]));
			}
			pass = 0;
		}
	}
	avg /= N*N*N;
	avgref /= N*N*N;
	printf("avg = (%g %g) ref = (%g %g)\n", real(avg), imag(avg), real(avgref), imag(avgref));
	if (pass)
		printf("PASSED!\n");
	else
		printf("FAILED! Total of %d errors found\n", errors);
	delete[] u;
	delete[] r;
	return 0;
}
