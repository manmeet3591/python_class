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

# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for similar SHT functions,
# from one generic function + tags.
# > See Makefile and SHT.c
# Basically, there are tags at the beginning of lines that are information
# to keep or remove the line depending on the function to build.
# tags :
# Q : line for scalar transform
# V : line for vector transform (both spheroidal and toroidal)
# S : line for vector transfrom, spheroidal component
# T : line for vector transform, toroidal component.

	static
3	void GEN3(_sy3,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, cplx *BrF, cplx *BtF, cplx *BpF, const long int llim, const int imlim)
QX	void GEN3(_sy1,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *BrF, const long int llim, const int imlim)
  #ifndef SHT_GRAD
VX	void GEN3(_sy2,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, cplx *BtF, cplx *BpF, const long int llim, const int imlim)
  #else
S	void GEN3(_sy1s,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *BtF, cplx *BpF, const long int llim, const int imlim)
T	void GEN3(_sy1t,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, cplx *BtF, cplx *BpF, const long int llim, const int imlim)
  #endif
  {
  #ifndef SHT_AXISYM
Q	#define qr(l) vall(creal(Ql[l-1]))
Q	#define qi(l) vall(cimag(Ql[l-1]))
S	#define sr(l) vall(creal(Sl[l-1]))
S	#define si(l) vall(cimag(Sl[l-1]))
T	#define tr(l) vall(creal(Tl[l-1]))
T	#define ti(l) vall(cimag(Tl[l-1]))
V	double m_1;
	unsigned im;
  #endif
	unsigned m0, mstep;
	long int nk,k,l,m;
	double *alm, *al;
	double *ct, *st;
Q	cplx Ql[llim+1] SSE;
S	cplx Sl[llim] SSE;
T	cplx Tl[llim] SSE;

Q	double rnr[NLAT_2 + NWAY*VSIZE2 -1] SSE;
Q	double rsr[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double tnr[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double tsr[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double pnr[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double psr[NLAT_2 + NWAY*VSIZE2 -1] SSE;
  #ifndef SHT_AXISYM
Q	double rni[NLAT_2 + NWAY*VSIZE2 -1] SSE;
Q	double rsi[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double tni[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double tsi[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double pni[NLAT_2 + NWAY*VSIZE2 -1] SSE;
V	double psi[NLAT_2 + NWAY*VSIZE2 -1] SSE;
  #endif


	ct = shtns->ct;		st = shtns->st;
	nk = NLAT_2;
	nk = ((unsigned)(nk+VSIZE2-1)) / VSIZE2;

	// ACCESS PATTERN
	const int k_inc = 1;		const int m_inc = NLAT/2;

	#ifndef _OPENMP
		m0 = 0;		mstep = 1;
	#else
		m0 = omp_get_thread_num();
		mstep = omp_get_num_threads();
		if (m0 == 0)
	#endif
	{	//	im=0;
		#ifdef SHT_GRAD
		  #ifndef SHT_AXISYM
S			k=0; do { BpF[k]=0.0; } while(++k<NLAT_2);
T			k=0; do { BtF[k]=0.0; } while(++k<NLAT_2);
		  #else
S			if (BpF != NULL) { int k=0; do { BpF[k]=0.0; } while(++k<NLAT_2); }
T			if (BtF != NULL) { int k=0; do { BtF[k]=0.0; } while(++k<NLAT_2); }
		  #endif
		#endif
Q		double* Ql0 = (double*) Ql;
S		double* Sl0 = (double*) Sl;
T		double* Tl0 = (double*) Tl;
 		l=1;
		alm = shtns->alm;
Q		Ql0[0] = (double) Qlm[0];		// l=0
		do {		// for m=0, compress the complex Q,S,T to double
Q			Ql0[l] = creal( Qlm[l] );	//	Ql[l+1] = (double) Qlm[l+1];
S			Sl0[l-1] = creal( Slm[l] );	//	Sl[l] = (double) Slm[l+1];
T			Tl0[l-1] = creal( Tlm[l] );	//	Tl[l] = (double) Tlm[l+1];
			++l;
		} while(l<=llim);
		k=0;
		do {
			l=0;	al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
V			rnd sint[NWAY], dy0[NWAY], dy1[NWAY];
Q			rnd re[NWAY], ro[NWAY];
S			rnd te[NWAY], to[NWAY];
T			rnd pe[NWAY], po[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, j+k);
V				sint[j] = -vread(st, j+k);
				y0[j] = vall(al[0]);
V				dy0[j] = vall(0.0);
Q				re[j] = y0[j] * vall(Ql0[0]);
S				to[j] = dy0[j];
T				po[j] = dy0[j];
			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = vall(al[0]*al[1]) * cost[j];
V				dy1[j] = vall(al[0]*al[1]) * sint[j];
			}
			for (int j=0; j<NWAY; ++j) {
Q				ro[j] = y1[j] * vall(Ql0[1]);
S				te[j] = dy1[j] * vall(Sl0[0]);
T				pe[j] = -dy1[j] * vall(Tl0[0]);
			}
			al+=2;	l+=2;
			while(l<llim) {
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					re[j] += y0[j] * vall(Ql0[l]);
S					to[j] += dy0[j] * vall(Sl0[l-1]);
T					po[j] -= dy0[j] * vall(Tl0[l-1]);
				}
				for (int j=0; j<NWAY; ++j) {
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*sint[j]) + vall(al[2])*dy1[j];
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					ro[j] += y1[j] * vall(Ql0[l+1]);
S					te[j] += dy1[j] * vall(Sl0[l]);
T					pe[j] -= dy1[j] * vall(Tl0[l]);
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*cost[j]*y1[j] + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					re[j] += y0[j] * vall(Ql0[l]);
S					to[j] += dy0[j] * vall(Sl0[l-1]);
T					po[j] -= dy0[j] * vall(Tl0[l-1]);
				}
			}

			for (int j=0; j<NWAY; ++j) {
Q				vstor(rnr, j+k, re[j]+ro[j]);		vstor(rsr, j+k, re[j]-ro[j]);
S				vstor(tnr, j+k, te[j]+to[j]);		vstor(tsr, j+k, te[j]-to[j]);
T				vstor(pnr, j+k, pe[j]+po[j]);		vstor(psr, j+k, pe[j]-po[j]);
			}
			k+=NWAY;
		} while (k < nk);

		k=0;  do {	// merge symmetric and antisymmetric parts.
Q			BrF[(k/2)*k_inc]            = rnr[k]   + I*rnr[k+1];
Q			BrF[(NLAT_2-1 - k/2)*k_inc] = rsr[k+1] + I*rsr[k];
S			BtF[(k/2)*k_inc]            = tnr[k]   + I*tnr[k+1];
S			BtF[(NLAT_2-1 - k/2)*k_inc] = tsr[k+1] + I*tsr[k];
T			BpF[(k/2)*k_inc]            = pnr[k]   + I*pnr[k+1];
T			BpF[(NLAT_2-1 - k/2)*k_inc] = psr[k+1] + I*psr[k];
			k+=2;
		} while(k < NLAT_2);

		m0=mstep;
	}

  #ifndef SHT_AXISYM
	for (im=m0; im<imlim; im+=mstep) {
		m = im*MRES;
V		m_1 = 1.0/m;
		//alm = shtns->alm[im];
		alm = shtns->alm + im*(2*LMAX -m+MRES);
		l = m;
		k = LiM(shtns, l,im);
		//k = (im*(2*(LMAX+1)-(m+MRES)))>>1 + l;
		do {		// copy input coefficients to a local array.
Q			((v2d*)Ql)[l-1] = ((v2d*)Qlm)[k];
S			((v2d*)Sl)[l-1] = ((v2d*)Slm)[k];
T			((v2d*)Tl)[l-1] = ((v2d*)Tlm)[k];
			++l;	++k;
		} while(l<=llim);

		k = shtns->tm[im] / VSIZE2;			// stay on vector boundary
		#if VSIZE2 == 1
			k -= k&1;		// we operate without vectors, but we still need complex alignement (2 doubles).
		#endif
		do {
			al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
V			rnd st2[NWAY], dy0[NWAY], dy1[NWAY];
Q			rnd rer[NWAY], rei[NWAY], ror[NWAY], roi[NWAY];
V			rnd ter[NWAY], tei[NWAY], tor[NWAY], toi[NWAY];
V			rnd per[NWAY], pei[NWAY], por[NWAY], poi[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(st, k+j);
				y0[j] = vall(1.0);
V				st2[j] = cost[j]*cost[j]*vall(-m_1);
V				y0[j] = vall(m);		// for the vector transform, compute ylm*m/sint
			}
Q			l=m;
V			l=m-1;
			long int ny = 0;
		if ((int)llim <= SHT_L_RESCALE_FLY) {
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
			} while(l >>= 1);
		} else {
			long int nsint = 0;
			do {		// sin(theta)^m		(use rescaling to avoid underflow)
				if (l&1) {
					for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
					ny += nsint;
					if (vlo(y0[0]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
						ny--;
						for (int j=0; j<NWAY; ++j) y0[j] *= vall(SHT_SCALE_FACTOR);
					}
				}
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
				nsint += nsint;
				if (vlo(cost[0]) < 1.0/SHT_SCALE_FACTOR) {
					nsint--;
					for (int j=0; j<NWAY; ++j) cost[j] *= vall(SHT_SCALE_FACTOR);
				}
			} while(l >>= 1);
		}
			for (int j=0; j<NWAY; ++j) {
				y0[j] *= vall(al[0]);
				cost[j] = vread(ct, j+k);
V				dy0[j] = cost[j]*y0[j];
Q				ror[j] = vall(0.0);		roi[j] = vall(0.0);
Q				rer[j] = vall(0.0);		rei[j] = vall(0.0);
			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = (vall(al[1])*y0[j]) *cost[j];		//	y1[j] = vall(al[1])*cost[j]*y0[j];
V				por[j] = vall(0.0);		tei[j] = vall(0.0);
V				tor[j] = vall(0.0);		pei[j] = vall(0.0);
V				dy1[j] = (vall(al[1])*y0[j]) *(cost[j]*cost[j] + st2[j]);		//	dy1[j] = vall(al[1])*(cost[j]*dy0[j] - y0[j]*st2[j]);
V				poi[j] = vall(0.0);		ter[j] = vall(0.0);
V				toi[j] = vall(0.0);		per[j] = vall(0.0);
			}
			l=m;		al+=2;
			while ((ny<0) && (l<llim)) {		// ylm treated as zero and ignored if ny < 0
				for (int j=0; j<NWAY; ++j) {
					y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*st2[j]) + vall(al[0])*dy0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*st2[j]) + vall(al[2])*dy1[j];				
				}
				l+=2;	al+=4;
				if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
					++ny;
					for (int j=0; j<NWAY; ++j) {
						y0[j] *= vall(1.0/SHT_SCALE_FACTOR);		y1[j] *= vall(1.0/SHT_SCALE_FACTOR);
V						dy0[j] *= vall(1.0/SHT_SCALE_FACTOR);		dy1[j] *= vall(1.0/SHT_SCALE_FACTOR);
					}
				}
			}
		  if (ny == 0) {
			while (l<llim) {	// compute even and odd parts
Q				for (int j=0; j<NWAY; ++j) {	rer[j] += y0[j]  * qr(l);		rei[j] += y0[j] * qi(l);	}
Q				for (int j=0; j<NWAY; ++j) {	ror[j] += y1[j]  * qr(l+1);		roi[j] += y1[j] * qi(l+1);	}
S				for (int j=0; j<NWAY; ++j) {	tor[j] += dy0[j] * sr(l);		pei[j] += y0[j] * sr(l);	}
S				for (int j=0; j<NWAY; ++j) {	ter[j] += dy1[j] * sr(l+1);		poi[j] += y1[j] * sr(l+1);	}
S				for (int j=0; j<NWAY; ++j) {	toi[j] += dy0[j] * si(l);		per[j] -= y0[j] * si(l);	}
S				for (int j=0; j<NWAY; ++j) {	tei[j] += dy1[j] * si(l+1);		por[j] -= y1[j] * si(l+1);	}
T				for (int j=0; j<NWAY; ++j) {	por[j] -= dy0[j] * tr(l);		tei[j] += y0[j] * tr(l);	}
T				for (int j=0; j<NWAY; ++j) {	per[j] -= dy1[j] * tr(l+1);		toi[j] += y1[j] * tr(l+1);	}
T				for (int j=0; j<NWAY; ++j) {	poi[j] -= dy0[j] * ti(l);		ter[j] -= y0[j] * ti(l);	}
T				for (int j=0; j<NWAY; ++j) {	pei[j] -= dy1[j] * ti(l+1);		tor[j] -= y1[j] * ti(l+1);	}
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*st2[j]) + vall(al[0])*dy0[j];
					y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*st2[j]) + vall(al[2])*dy1[j];
					y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				l+=2;	al+=4;
			}
			if (l==llim) {
Q				for (int j=0; j<NWAY; ++j) {	rer[j] += y0[j] * qr(l);		rei[j] += y0[j] * qi(l);	}
S				for (int j=0; j<NWAY; ++j) {	tor[j] += dy0[j] * sr(l);		pei[j] += y0[j] * sr(l);	}
S				for (int j=0; j<NWAY; ++j) {	toi[j] += dy0[j] * si(l);		per[j] -= y0[j] * si(l);	}
T				for (int j=0; j<NWAY; ++j) {	por[j] -= dy0[j] * tr(l);		tei[j] += y0[j] * tr(l);	}
T				for (int j=0; j<NWAY; ++j) {	poi[j] -= dy0[j] * ti(l);		ter[j] -= y0[j] * ti(l);	}
			}
3			for (int j=0; j<NWAY; ++j) cost[j]  = vread(st, k+j) * vall(m_1);
3			for (int j=0; j<NWAY; ++j) {  rer[j] *= cost[j];  ror[j] *= cost[j];	rei[j] *= cost[j];  roi[j] *= cost[j];  }
		  }
		  
			for (int j=0; j<NWAY; ++j) {
Q				vstor(rnr, j+k, rer[j]+ror[j]);		vstor(rsr, j+k, rer[j]-ror[j]);
Q				vstor(rni, j+k, rei[j]+roi[j]);		vstor(rsi, j+k, rei[j]-roi[j]);
V				vstor(tnr, j+k, ter[j]+tor[j]);		vstor(tsr, j+k, ter[j]-tor[j]);
V				vstor(tni, j+k, tei[j]+toi[j]);		vstor(tsi, j+k, tei[j]-toi[j]);
V				vstor(pnr, j+k, per[j]+por[j]);		vstor(psr, j+k, per[j]-por[j]);
V				vstor(pni, j+k, pei[j]+poi[j]);		vstor(psi, j+k, pei[j]-poi[j]);
			}
			k+=NWAY;
		} while (k < nk);

		l = shtns->tm[im] >> 1;		// stay on a 16 byte boundary
Q		k=0;	while (k<l) {	// polar optimization
Q			BrF[im*m_inc + k*k_inc] = 0.0;				BrF[(NPHI-im)*m_inc + k*k_inc] = 0.0;
Q			BrF[im*m_inc + (NLAT_2-l+k)*k_inc] = 0.0;	BrF[(NPHI-im)*m_inc + (NLAT_2-l+k)*k_inc] = 0.0;
Q			++k;
Q		}
Q		k*=2;	do {
Q			BrF[im*m_inc + (k/2)*k_inc] = (rnr[k]-rni[k+1]) + I*(rnr[k+1]+rni[k]);
Q			BrF[(NPHI-im)*m_inc + (k/2)*k_inc] = (rnr[k]+rni[k+1]) + I*(rnr[k+1]-rni[k]);
Q			BrF[im*m_inc + (NLAT_2-1-k/2)*k_inc] = (rsr[k+1]-rsi[k]) + I*(rsr[k]+rsi[k+1]);
Q			BrF[(NPHI-im)*m_inc + (NLAT_2-1-k/2)*k_inc] = (rsr[k+1]+rsi[k]) + I*(rsr[k]-rsi[k+1]);
Q			k+=2;
Q		} while(k < NLAT_2);

V		k=0;	while (k<l) {	// polar optimization
V			BtF[im*m_inc + k*k_inc] = 0.0;				BtF[(NPHI-im)*m_inc + k*k_inc] = 0.0;
V			BtF[im*m_inc + (NLAT_2-l+k)*k_inc] = 0.0;	BtF[(NPHI-im)*m_inc + (NLAT_2-l+k)*k_inc] = 0.0;
V			++k;
V		}
V		k*=2;	do {
V			BtF[im*m_inc + (k/2)*k_inc] = (tnr[k]-tni[k+1]) + I*(tnr[k+1]+tni[k]);
V			BtF[(NPHI-im)*m_inc + (k/2)*k_inc] = (tnr[k]+tni[k+1]) + I*(tnr[k+1]-tni[k]);
V			BtF[im*m_inc + (NLAT_2-1-k/2)*k_inc] = (tsr[k+1]-tsi[k]) + I*(tsr[k]+tsi[k+1]);
V			BtF[(NPHI-im)*m_inc + (NLAT_2-1-k/2)*k_inc] = (tsr[k+1]+tsi[k]) + I*(tsr[k]-tsi[k+1]);
V			k+=2;
V		} while(k < NLAT_2);

V		k=0;	while (k<l) {	// polar optimization
V			BpF[im*m_inc + k*k_inc] = 0.0;				BpF[(NPHI-im)*m_inc + k*k_inc] = 0.0;
V			BpF[im*m_inc + (NLAT_2-l+k)*k_inc] = 0.0;	BpF[(NPHI-im)*m_inc + (NLAT_2-l+k)*k_inc] = 0.0;
V			++k;
V		}
V		k*=2;	do {
V			BpF[im*m_inc + (k/2)*k_inc] = (pnr[k]-pni[k+1]) + I*(pnr[k+1]+pni[k]);
V			BpF[(NPHI-im)*m_inc + (k/2)*k_inc] = (pnr[k]+pni[k+1]) + I*(pnr[k+1]-pni[k]);
V			BpF[im*m_inc + (NLAT_2-1-k/2)*k_inc] = (psr[k+1]-psi[k]) + I*(psr[k]+psi[k+1]);
V			BpF[(NPHI-im)*m_inc + (NLAT_2-1-k/2)*k_inc] = (psr[k+1]+psi[k]) + I*(psr[k]-psi[k+1]);
V			k+=2;
V		} while(k < NLAT_2);

	}

	while(im <= NPHI-imlim) {	// padding for high m's
		k=0;
		do {
Q			BrF[im*m_inc + k*k_inc] = 0.0;
V			BtF[im*m_inc + k*k_inc] = 0.0;
V			BpF[im*m_inc + k*k_inc] = 0.0;
		} while (++k < NLAT_2);
	  im+=mstep;
	}
  #endif
  }

Q	#undef qr
Q	#undef qi
S	#undef sr
S	#undef si
T	#undef tr
T	#undef ti

	static
3	void GEN3(SHqst_to_spat_mic,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, long int llim) {
QX	void GEN3(SH_to_spat_mic,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, double *Vr, long int llim) {
  #ifndef SHT_GRAD
VX	void GEN3(SHsphtor_to_spat_mic,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, long int llim) {
  #else
S	void GEN3(SHsph_to_spat_mic,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp, long int llim) {
T	void GEN3(SHtor_to_spat_mic,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp, long int llim) {
  #endif

	int k;
	unsigned imlim = 0;
Q	cplx* BrF = (cplx*) Vr;
V	cplx* BtF = (cplx*) Vt;	cplx* BpF = (cplx*) Vp;

  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	if (shtns->fftc_mode > 0) {		// alloc memory for the FFT
		unsigned long nv = shtns->nspat;
QX		BrF = (cplx*) VMALLOC( nv * sizeof(double) );
VX		BtF = (cplx*) VMALLOC( 2*nv * sizeof(double) );
VX		BpF = BtF + nv/2;
3		BrF = (cplx*) VMALLOC( 3*nv * sizeof(double) );
3		BtF = BrF + nv/2;		BpF = BrF + nv;
	}
  #endif
	imlim += 1;
  
  #pragma omp parallel num_threads(shtns->nthreads)
  {
3	GEN3(_sy3,NWAY,SUFFIX)(shtns, Qlm, Slm, Tlm, BrF, BtF, BpF, llim, imlim);
QX	GEN3(_sy1,NWAY,SUFFIX)(shtns, Qlm, BrF, llim, imlim);
	#ifndef SHT_GRAD
VX		GEN3(_sy2,NWAY,SUFFIX)(shtns, Slm, Tlm, BtF, BpF, llim, imlim);
	#else
S		GEN3(_sy1s,NWAY,SUFFIX)(shtns, Slm, BtF, BpF, llim, imlim);
T		GEN3(_sy1t,NWAY,SUFFIX)(shtns, Tlm, BtF, BpF, llim, imlim);
	#endif
  }

  #ifndef SHT_AXISYM
    // NPHI > 1 as SHT_AXISYM is not defined.
  	if (shtns->fftc_mode >= 0) {
		if (shtns->fftc_mode == 0) {
Q			fftw_execute_dft(shtns->ifftc, (cplx *) BrF, (cplx *) Vr);
V			fftw_execute_dft(shtns->ifftc, (cplx *) BtF, (cplx *) Vt);
V			fftw_execute_dft(shtns->ifftc, (cplx *) BpF, (cplx *) Vp);
		} else {		// split dft
Q			fftw_execute_split_dft(shtns->ifftc,((double*)BrF)+1, ((double*)BrF), Vr+NPHI, Vr);
V			fftw_execute_split_dft(shtns->ifftc,((double*)BtF)+1, ((double*)BtF), Vt+NPHI, Vt);
V			fftw_execute_split_dft(shtns->ifftc,((double*)BpF)+1, ((double*)BpF), Vp+NPHI, Vp);
Q			VFREE(BrF);
VX			VFREE(BtF);		// this frees also BpF.
		}
	}
  #endif

  }
