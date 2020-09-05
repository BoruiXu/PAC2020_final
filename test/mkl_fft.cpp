#include "mkl_dfti.h"

float _Complex c2c_data[32];
DFTI_DESCRIPTOR_HANDLE my_desc1_handle = NULL;
MKL_LONG status;



/* ...put values into c2c_data[i] 0<=i<=31 */
/* ...put values into r2c_data[i] 0<=i<=31 */

for(int i=0;i<32;i++){
	c2c_data[i] = complex<float>(1,1);
}

status = DftiCreateDescriptor(&my_desc1_handle, DFTI_SINGLE,
                              DFTI_COMPLEX, 1, 32);
status = DftiCommitDescriptor(my_desc1_handle);
status = DftiComputeForward(my_desc1_handle, c2c_data);
status = DftiFreeDescriptor(&my_desc1_handle);