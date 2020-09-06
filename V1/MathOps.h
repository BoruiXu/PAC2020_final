#ifndef MATHSTUFF_H
#define MATHSTUFF_H

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

/* Create a L-vector of type T between a and b (inclusive) */
template<typename T> T* linspace(T a,T b,int L) {
	T *f = new T[L];
	for (int l=0;l<L;l++) f[l] = a+(b-a)*((T)l/(L-1));
	return f;
}

/* Modulus function (Knuth-style) */
template<typename T> T mod(T x,T n) {
	return x - n*(T)floor(1.0*x/n);
}

/* Round function (symmetric) */
template<typename T> T round(T x) {
	T temp;
	if (x>=0) {
		temp = x-floor(x);
		if (temp>=0.5) return ceil(x);
		else return floor(x);
	} else {
		temp = x-ceil(x);
		if (temp>=-0.5) return ceil(x);
		else return floor(x);
	}
}

/* 3D Chop */
template<typename T> void chop3D(T* f,int Nx,int Ny,int Nz,int chopX,int chopY,int chopZ) {
	if (chopX==1) {
		for (int x=0;x<Nx;x+=2) {
			for (int y=0;y<Ny;y++) {
				for (int z=0;z<Nz;z++) {
					f[x*Ny*Nz+y*Nz+z] = -f[x*Ny*Nz+y*Nz+z];
				}
			}
		}	
	}
	if (chopY==1) {
		for (int x=0;x<Nx;x++) {
			for (int y=0;y<Ny;y+=2) {
				for (int z=0;z<Nz;z++) {
					f[x*Ny*Nz+y*Nz+z] = -f[x*Ny*Nz+y*Nz+z];
				}
			}
		}
	}
	if (chopZ==1) {
		for (int x=0;x<Nx;x++) {
			for (int y=0;y<Ny;y++) {
				for (int z=0;z<Nz;z+=2) {
					f[x*Ny*Nz+y*Nz+z] = -f[x*Ny*Nz+y*Nz+z];	
				}
			}
		}
	}
	return;
}

#endif