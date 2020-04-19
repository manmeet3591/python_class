/*
 * Copyright (c) 2010-2015 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@ujf-grenoble.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */


/// \example test_rot.c
/// This program test the rotation of spherical harmonics (beta)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#include <shtns.h>

shtns_cfg shtns;

complex double *Slm, *Slm0, *Tlm, *Tlm0, *Qlm;	// spherical harmonics l,m space
complex double *ShF, *ThF, *NLF;	// Fourier space : theta,m
double *Sh, *Th, *NL;		// real space : theta,phi (alias of ShF)

int LMAX,MMAX,MRES,NLM;
int NLAT = 0;
int NPHI = 0;

void write_shell(char *fn, double *V)
{
	FILE *fp;
	int i,j,k;	//phi, theta
	
	fp = fopen(fn,"w");
	fprintf(fp,"0 ");
		for(j=0;j<NLAT;j++) {
			fprintf(fp,"%.6g 0 0 ",acos(shtns->ct[j]));	// first line = theta (radians)
		}
	for (i=0; i<NPHI; i++) {
		fprintf(fp,"\n%.6g ",PHI_RAD(shtns, i));		// first row = phi (radians)
		for(j=0; j<NLAT; j++) {
			k = i*NLAT+j;
			fprintf(fp,"%.6g %.6g %.6g  ",V[k],0.,0.);		// data
		}
	}
	fprintf(fp,"\n");	fclose(fp);
}

int main(int argc, char *argv[])
{
	complex double t1, t2;
	double t,tmax,n2;
	int i,im,m,l,jj;
	double e0,e1;
	double polaropt = 1.e-8;		// default for polar optimization.
	enum shtns_type shtmode = sht_auto;		// default to "auto" (fastest) mode.
	enum shtns_norm shtnorm = sht_orthonormal;		// default to "orthonormal" SH.
	int layout = SHT_NATIVE_LAYOUT;
	int nlorder = 0;
	int vector = 1;
	char name[20];

	srand( time(NULL) );	// initialise les nombres.
	
	MMAX=LMAX=7;
	MRES=1;
	NLAT=12;
	NPHI=16;

	shtns_use_threads(0);
	shtns = shtns_create(LMAX, MMAX, MRES, shtnorm);
	NLM = shtns->nlm;
	shtns_set_grid_auto(shtns, sht_quick_init, 0.0, 1e-10, &NLAT, &NPHI);
	shtns_print_cfg(shtns);

	ShF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Sh = (double *) ShF;
	ThF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	Th = (double *) ThF;
	NLF = (complex double *) fftw_malloc( 4*(NPHI/2+1) * NLAT * sizeof(complex double));
	NL = (double *) NLF;

	Tlm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm0 = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Slm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Tlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);
	Qlm = (complex double *) fftw_malloc(sizeof(complex double)* NLM);


	double mx[3][3];

	for (l=0; l<NLM; l++) 	Qlm[l] = 0.0;
	Qlm[LiM(shtns, 1, 0)] = 1.0;
	SH_Xrotate90(shtns, Qlm, Slm);
	mx[0][0] = Slm[LiM(shtns, 1, 0)];
	mx[1][0] = creal(Slm[LiM(shtns, 1, 1)]);
	mx[2][0] = cimag(Slm[LiM(shtns, 1, 1)]);

	for (l=0; l<NLM; l++) 	Qlm[l] = 0.0;
	Qlm[LiM(shtns, 1, 1)] = 1.0;
	SH_Xrotate90(shtns, Qlm, Slm);
	mx[0][1] = Slm[LiM(shtns, 1, 0)];
	mx[1][1] = creal(Slm[LiM(shtns, 1, 1)]);
	mx[2][1] = cimag(Slm[LiM(shtns, 1, 1)]);

	for (l=0; l<NLM; l++) 	Qlm[l] = 0.0;
	Qlm[LiM(shtns, 1, 1)] = I;
	SH_Xrotate90(shtns, Qlm, Slm);
	mx[0][2] = Slm[LiM(shtns, 1, 0)];
	mx[1][2] = creal(Slm[LiM(shtns, 1, 1)]);
	mx[2][2] = cimag(Slm[LiM(shtns, 1, 1)]);

	printf("rotation matrix Xrotate90:\n");
	printf(" %10f %10f %10f\n", mx[0][0], mx[0][1], mx[0][2]);
	printf(" %10f %10f %10f\n", mx[1][0], mx[1][1], mx[1][2]);
	printf(" %10f %10f %10f\n", mx[2][0], mx[2][1], mx[2][2]);
	printf("\n");
	
	for (l=0; l<NLM; l++) 	{Qlm[l] = 0.0;	Slm[l] = 0.0;}
//	Qlm[0] = 1000.;
//	Qlm[LiM(shtns, 2, 0)] = 1.0;
	Qlm[LiM(shtns, 1, 1)] = I;	//0.1 + I*0.05;
	
	SH_to_spat(shtns, Qlm, Sh);
	write_shell("q0",Sh);

	SH_Yrotate(shtns, Qlm, 30*M_PI/180., Slm);
	SH_to_spat(shtns, Slm, Sh);
	write_shell("qry30",Sh);

	SH_Yrotate90(shtns, Qlm, Slm);
	SH_to_spat(shtns, Slm, Sh);
	write_shell("qry",Sh);

	SH_Xrotate90(shtns, Qlm, Slm);
	SH_to_spat(shtns, Slm, Sh);
	write_shell("qrx",Sh);

	SH_Zrotate(shtns, Qlm, 30*M_PI/180, Slm);
	SH_to_spat(shtns, Slm, Sh);
	write_shell("q0z",Sh);


	printf("\n*** rotation along Z by 30 deg.\n");
	for (l=0; l<=LMAX; l++) {
		double norm0=0;		double normR=0;
		for (m=0; m<=l; m++) {
			printf("l=%d, m=%d,  Q=%f,%f,  \t  S=%f,%f\n",l,m, creal(Qlm[LiM(shtns,l,m)]), cimag(Qlm[LiM(shtns,l,m)]),
			creal(Slm[LiM(shtns,l,m)]), cimag(Slm[LiM(shtns,l,m)]));
			double e0 = creal(Qlm[LiM(shtns,l,m)])*creal(Qlm[LiM(shtns,l,m)]) + cimag(Qlm[LiM(shtns,l,m)])*cimag(Qlm[LiM(shtns,l,m)]);
			double eR = creal(Slm[LiM(shtns,l,m)])*creal(Slm[LiM(shtns,l,m)]) + cimag(Slm[LiM(shtns,l,m)])*cimag(Slm[LiM(shtns,l,m)]);			
			if (m>0) {
				e0 *=2;		eR *= 2;
			}
			norm0 += e0;
			normR += eR;
		}
		printf("norm0 = %f, norm1 = %f\n", sqrt(norm0), sqrt(normR));
	}


//	SH_Yrotate90(shtns, Slm, Tlm0);
//	SH_Yrotate90(shtns, Tlm0, Slm);
	SH_Yrotate90(shtns, Slm, Tlm);

	printf("\n*** rotation along Y by 90 deg.\n");
	for (l=0; l<=LMAX; l++) {
		double norm0=0;		double normR=0;
		for (m=0; m<=l; m++) {
			printf("l=%d, m=%d,  Q=%f,%f,  \t  S=%f,%f\n",l,m, creal(Qlm[LiM(shtns,l,m)]), cimag(Qlm[LiM(shtns,l,m)]),
			creal(Tlm[LiM(shtns,l,m)]), cimag(Tlm[LiM(shtns,l,m)]));
			double e0 = creal(Qlm[LiM(shtns,l,m)])*creal(Qlm[LiM(shtns,l,m)]) + cimag(Qlm[LiM(shtns,l,m)])*cimag(Qlm[LiM(shtns,l,m)]);
			double eR = creal(Tlm[LiM(shtns,l,m)])*creal(Tlm[LiM(shtns,l,m)]) + cimag(Tlm[LiM(shtns,l,m)])*cimag(Tlm[LiM(shtns,l,m)]);			
			if (m>0) {
				e0 *=2;		eR *= 2;
			}
			norm0 += e0;
			normR += eR;
		}
		printf("norm0 = %f, norm1 = %f\n", sqrt(norm0), sqrt(normR));
	}

}

