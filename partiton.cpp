
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <x86intrin.h>
#include <random>

using namespace std;

int hist(int end , float* array,int len){
	int i=0;
	while(array[i]<=end){
		i++;
	}
	return i;
}

void SORT(float* x,float* y,float* z,int d,int len){
	float* priority;
	float* t1;
	float* t2;
	if(d == 0){
		priority = x;
		t1 = y;
		t2 = z;
	}
	else if(d == 1){
		priority = y;
		t1 = x;
		t2 = z;
	}
	else if(d == 2){
		priority = z;
		t1 = y;
		t2 = x;
	}

	for(int i = 0; i<len;i++){
		for(int j=i+1;j<len;j++){
			if(priority[i]>priority[j]){
				float temp = priority[i];
				priority[i] = priority[j];
				priority[j] = temp;

				temp = t1[i];
				t1[i] = t1[j];
				t1[j] = temp;

				temp = t2[i];
				t2[i] = t2[j];
				t2[j] = temp;
			}
		}
	}

}
void task(int M,float* wx, float* wy,float* wz,int P, int W,float** part,int* numPart){
	int average = M/P;
	int minWidth = 2*W+1;
	//int part[3][P];
	//int numPart[3];
	for(int i=0;i<3;i++){
		//sort
		SORT(wx,wy,wz,i,M);
		// for(int j = 0;j<M;j++){
		// 	if(i==0){
		// 		cout<<wx[j]<<" ";
		// 	}
		// 	else if(i==1){
		// 		cout<<wy[j]<<" ";
		// 	}
		// 	else{
		// 		cout<<wz[j]<<" ";
		// 	}
		// }
		// cout<<endl;

		int index =0;
		float end,start,N;
		float* array;
		if(i==0){
			array = wx;
		}
		else if(i==1){
			array = wy;
		}
		else{
			array = wz;
		}


		start = array[0];
		N = array[M-1];
		
		part[i][0] = start;
		while(start<N){

			end = start+minWidth;
			if(end>N){
				end = N;
			}
			
			
			while(hist(end,array,M)-hist(start,array,M)<average){
				end++;
				if(end>=N){
					end=N;
					break;
				}
			}

			part[i][++index] = end;
			start = end;

		}

		numPart[i] = index;

	}
}




int main(){
	float x[10] = {100.001,144.581,125.749,137.176,129.111,125.586,136.311,121.303,145.077,142.899};
	float y[10] = {104.252,148.398,119.9,104.477,140.478,143.832,148.331,144.975,148.077,145.342};
	float z[10] = {130.068,109.484,113.145,128.02,129.596,149.754,114.855,132.65,108.236,114.701};
	// default_random_engine e;
 //    uniform_real_distribution<float> u(100, 150);

	// for(int i=0;i<10;i++){
	// 	x[i] = u(e);
	// 	y[i] = u(e);
	// 	z[i] = u(e);
	// }
	// for(int i= 0 ;i<10;i++){
	// 	cout<<x[i]<<",";
	// }
	// cout<<endl;
	
	// for(int i= 0 ;i<10;i++){
	// 	cout<<y[i]<<",";
	// }
	// cout<<endl;
	
	// for(int i= 0 ;i<10;i++){
	// 	cout<<z[i]<<",";
	// }
	// cout<<endl;




	float** part;
	part = (float**)malloc(sizeof(float*)*(3));
	for(int i=0;i<3;i++){
		part[i] = (float*)malloc((10)*sizeof(float));
	}  
	int numPart[3];
	task(10,x,y,z,5,2,part,numPart);
	for(int i=0;i<3;i++){
		cout<<numPart[i]<<endl;
	}
	for(int j =0;j<3;j++){
		for(int i=0;i<6;i++){
			cout<<part[j][i]<<" ";
		}
		cout<<endl;
	}
	

	

	for(int i=0;i<3;i++){
		 free((void *)part[i]);
	}
	free((void *)part);





}
























