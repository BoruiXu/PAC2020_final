#include "mkl_dfti.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <complex.h>
#include <malloc.h>

using namespace std;

int main(){
	//float _Complex c2c_data[32];


	complex<float> *c2c_data;
	c2c_data=(complex<float>*)memalign(16,32*sizeof(complex<float>));
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle = NULL;
	MKL_LONG status;



	/* ...put values into c2c_data[i] 0<=i<=31 */
	/* ...put values into r2c_data[i] 0<=i<=31 */

	for(int i=0;i<32;i++){
		c2c_data[i] = complex<float>(2,5);
	}
	

	status = DftiCreateDescriptor(&my_desc1_handle, DFTI_SINGLE,
	                              DFTI_COMPLEX, 1, 32);
	status = DftiCommitDescriptor(my_desc1_handle);
	status = DftiComputeForward(my_desc1_handle, c2c_data);
	status = DftiFreeDescriptor(&my_desc1_handle);
	for(int i=0;i<32;i++){
		//printf("%f+%fi\n",creal(c2c_data[i]),cimag(c2c_data[i]));
		std::cout << "Real part: " << c2c_data[i].real()<< " Imaginary part: " << c2c_data[i].imag()<< std::endl;
	}


	return 0;
}
