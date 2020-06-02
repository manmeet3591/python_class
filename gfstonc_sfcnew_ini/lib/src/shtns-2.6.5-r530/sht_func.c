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

/** \internal \file sht_rot.c
 * \brief Rotation of Spherical Harmonics.
 */


/** \addtogroup rotation Rotation of SH fields.
Rotation around axis other than Z should be considered of beta quality (they have been tested but may still contain bugs).
They also require \c mmax = \c lmax. They use an Algorithm inspired by the pseudospectral rotation described in
Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" <i>Journal of Computational Physics</i>.
doi:<a href="http://dx.doi.org/10.1016/j.jcp.2009.05.014">10.1016/j.jcp.2009.05.014</a>

These functions do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
*/
//@{

/// Rotate a SH representation Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	int im, l, lmax, mmax, mres;

	lmax = shtns->lmax;		mmax = shtns->mmax;		mres = shtns->mres;

	if (Rlm != Qlm) {		// copy m=0 which does not change.
		l=0;	do { Rlm[l] = Qlm[l]; } while(++l <= lmax);
	}
	if (mmax > 0) {
		im=1; do {
			cplx eima = cos(im*mres*alpha) - I*sin(im*mres*alpha);		// rotate reference frame by angle -alpha
			for (l=im*mres; l<=lmax; ++l)	Rlm[LiM(shtns, l, im)] = Qlm[LiM(shtns, l, im)] * eima;
		} while(++im <= mmax);
	}
}

//@}

/** \internal rotation kernel used by SH_Yrotate90(), SH_Xrotate90() and SH_rotate().
 Algorithm based on the pseudospectral rotation[1] :
 - rotate around Z by angle dphi0.
 - synthetize for each l the spatial description for phi=0 and phi=pi on an equispaced latitudinal grid.
 - Fourier ananlyze as data on the equator to recover the m in the 90 degrees rotated frame.
 - rotate around new Z by angle dphi1.
 [1] Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" Journal of Computational Physics. **/
static void SH_rotK90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm, double dphi0, double dphi1)
{
	fftw_plan fft;
	cplx *q;
	double *q0;
	long int k, m, l;
	int lmax, ntheta;
	int nrembed, ncembed;

	lmax = shtns->lmax;
	ntheta = ((lmax+2)>>1)*2;
	m = 2* sizeof(double)*(2*ntheta+2)*lmax;
	q0 = VMALLOC(m);		memset(q0, 0, m);		// alloc & zero out.

	// rotate around Z by dphi0
	if (dphi0 != 0.0) {
		SH_Zrotate(shtns, Qlm, dphi0, Rlm);
		Qlm = Rlm;
	} else {
		Rlm[0] = Qlm[0];		// l=0 is rotation invariant.
	}

  #pragma omp parallel private(k,m,l) num_threads(shtns->nthreads)
  {
	double yl[lmax+1];
	// compute q(l) on the meridian phi=0 and phi=pi. (rotate around X)
	#pragma omp for schedule(static)
	for (k=0; k<ntheta/2; ++k) {
		double cost= cos(((0.5*M_PI)*(2*k+1))/ntheta);
		double sint_1 = 1.0/sqrt((1.0-cost)*(1.0+cost));
		m=0;
			legendre_sphPlm_array(shtns, lmax, m, cost, yl+m);
			double sgnt = -1.0;
			for (l=1; l<=lmax; ++l) {
				double qr = creal(Qlm[LiM(shtns, l, m)]) * yl[l];
				q0[k*2*lmax +2*(l-1)] = qr;
				q0[(ntheta-1-k)*2*lmax +2*(l-1)] = sgnt*qr;
				q0[(ntheta+k)*2*lmax +2*(l-1)] = sgnt*qr;
				q0[(2*ntheta-1-k)*2*lmax +2*(l-1)] = qr;
				sgnt *= -1.0;
			}
		#if _GCC_VEC_ && __SSE2__
		s2d sgnm = SIGN_MASK_HI;
		s2d sgnflip = SIGN_MASK_2;
		for (m=1; m<=lmax; ++m) {
			legendre_sphPlm_array(shtns, lmax, m, cost, yl+m);
			s2d sgnt = vdup(0.0);
			s2d m_st = vset(2.0, -2*m*sint_1);		// x2 for m>0
			sgnm = _mm_xor_pd(sgnm, sgnflip);	// (-1)^m
			for (l=m; l<=lmax; ++l) {
				v2d qc = ((v2d*)Qlm)[LiM(shtns, l, m)] * vdup(yl[l]) * m_st;	// (q0, dq0)
				((v2d*)q0)[k*lmax +(l-1)] += qc;
				((v2d*)q0)[(ntheta-1-k)*lmax +(l-1)] += (v2d)_mm_xor_pd(sgnt, qc);
				qc = _mm_xor_pd(sgnm, qc);
				((v2d*)q0)[(ntheta+k)*lmax +(l-1)] += (v2d)_mm_xor_pd( sgnt, qc );
				((v2d*)q0)[(2*ntheta-1-k)*lmax +(l-1)] += qc;
				sgnt = _mm_xor_pd(sgnt, sgnflip);	// (-1)^(l+m)
			}
		}
		#else
		double sgnm = 1.0;
		for (m=1; m<=lmax; ++m) {
			legendre_sphPlm_array(shtns, lmax, m, cost, yl+m);
			double sgnt = 1.0;
			sgnm *= -1.0;
			for (l=m; l<=lmax; ++l) {
				double qr = creal(Qlm[LiM(shtns, l, m)]) * yl[l];
				double qi = cimag(Qlm[LiM(shtns, l, m)]) * m*yl[l]*sint_1;
				qr += qr;	qi += qi;		// x2 for m>0
				q0[k*2*lmax +2*(l-1)] += qr;						// q0
				q0[k*2*lmax +2*(l-1)+1] -= qi;						// dq0
				q0[(ntheta-1-k)*2*lmax +2*(l-1)] += sgnt*qr;
				q0[(ntheta-1-k)*2*lmax +2*(l-1)+1] -= sgnt*qi;
				q0[(ntheta+k)*2*lmax +2*(l-1)] += (sgnm*sgnt)*qr;
				q0[(ntheta+k)*2*lmax +2*(l-1)+1] += (sgnm*sgnt)*qi;
				q0[(2*ntheta-1-k)*2*lmax +2*(l-1)] += sgnm*qr;
				q0[(2*ntheta-1-k)*2*lmax +2*(l-1)+1] += sgnm*qi;
				sgnt *= -1.0;
			}
		}
		#endif
	}
  }

	// perform FFT
	#ifdef OMP_FFTW
		k = (lmax < 63) ? 1 : shtns->nthreads;
		fftw_plan_with_nthreads(k);
	#endif
	q = (cplx*) q0;
	ntheta*=2;		nrembed = ntheta+2;		ncembed = nrembed/2;
	fft = fftw_plan_many_dft_r2c(1, &ntheta, 2*lmax, q0, &nrembed, 2*lmax, 1, q, &ncembed, 2*lmax, 1, FFTW_ESTIMATE);
	fftw_execute_dft_r2c(fft, q0, q);
	fftw_destroy_plan(fft);

	double yl[lmax+1];		double dyl[lmax+1];
	m=0;
		//legendre_sphPlm_deriv_array(shtns, lmax, m, 0.0, 1.0, yl+m, dyl+m);
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, yl+m, dyl+m);
		for (l=1; l<lmax; l+=2) {
			Rlm[LiM(shtns, l,m)] =  -creal(q[m*2*lmax +2*(l-1)+1])/(dyl[l]*ntheta);
			Rlm[LiM(shtns, l+1,m)] =  creal(q[m*2*lmax +2*l])/(yl[l+1]*ntheta);
		}
		if (l==lmax) {
			Rlm[LiM(shtns, l,m)] =  -creal(q[m*2*lmax +2*(l-1)+1])/(dyl[l]*ntheta);
		}
	dphi1 += M_PI/ntheta;	// shift rotation angle by angle of first synthesis latitude.
	for (m=1; m<=lmax; ++m) {
		//legendre_sphPlm_deriv_array(shtns, lmax, m, 0.0, 1.0, yl+m, dyl+m);
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, yl+m, dyl+m);
		cplx eimdp = (cos(m*dphi1) - I*sin(m*dphi1))/(ntheta);
		for (l=m; l<lmax; l+=2) {
			Rlm[LiM(shtns, l,m)] =  eimdp*q[m*2*lmax +2*(l-1)]*(1./yl[l]);
			Rlm[LiM(shtns, l+1,m)] =  eimdp*q[m*2*lmax +2*l+1]*(-1./dyl[l+1]);
		}
		if (l==lmax) {
			Rlm[LiM(shtns, l,m)] =  eimdp*q[m*2*lmax +2*(l-1)]*(1./yl[l]);
		}
	}
	VFREE(q0);
}


/// \addtogroup rotation
//@{

/// rotate Qlm by 90 degrees around X axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Xrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;													// rotation matrix for rotX(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * cimag(Qlm[LiM(shtns, l, 1)]);			//[m=0]     0        0    sqrt(2)
			Rlm[LiM(shtns, l ,1)] = creal(Qlm[LiM(shtns, l, 1)]) - I*(sqrt(0.5)*q0);	//[m=1r]    0        1      0
		return;																			//[m=1i] -sqrt(2)/2  0      0
	}

	SH_rotK90(shtns, Qlm, Rlm, 0.0,  -M_PI/2);
}

/// rotate Qlm by 90 degrees around Y axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Yrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;											// rotation matrix for rotY(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);									//[m=0]       0       0    sqrt(2)
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * creal(Qlm[LiM(shtns, l, 1)]);			//[m=1r] -sqrt(2)/2   0      0
			Rlm[LiM(shtns, l ,1)] = I*cimag(Qlm[LiM(shtns, l, 1)]) - sqrt(0.5) * q0;	//[m=1i]      0       0      1
		return;
	}

	SH_rotK90(shtns, Qlm, Rlm, -M_PI/2, 0.0);
}

/// rotate Qlm around Y axis by arbitrary angle, using composition of rotations. Store the result in Rlm.
void SH_Yrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if ((shtns->mres != 1) || (shtns->mmax < shtns->lmax)) shtns_runerr("truncature makes rotation not closed.");

	SH_rotK90(shtns, Qlm, Rlm, 0.0, M_PI/2 + alpha);	// Zrotate(pi/2) + Yrotate90 + Zrotate(pi+alpha)
	SH_rotK90(shtns, Rlm, Rlm, 0.0, M_PI/2);			// Yrotate90 + Zrotate(pi/2)
}

//@}



/** \addtogroup operators Special operators
 * Apply special operators in spectral space: multiplication by cos(theta), sin(theta).d/dtheta.
*/
//@{




/// fill mx with the coefficients for multiplication by cos(theta)
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// xq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void mul_ct_matrix(shtns_cfg shtns, double* mx)
{
	long int im,l,lm;
	double a_1;

	if (SHT_NORM == sht_schmidt) {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = alm_im(shtns,im);
			long int m=im*MRES;
			mx[2*lm] = 0.0;
			a_1 = 1.0 / al[1];
			l=m;
			while(++l < LMAX) {
				al+=2;				
				mx[2*lm+2] = a_1;
				a_1 = 1.0 / al[1];
				mx[2*lm+1] = -a_1*al[0];        // = -al[2*(lm+1)] / al[2*(lm+1)+1];
				lm++;
			}
			if (l == LMAX) {	// the last one needs to be computed.
				mx[2*lm+2] = a_1;
				mx[2*lm+1] = sqrt((l+m)*(l-m))/(2*l+1);
				lm++;
			}
			mx[2*lm +1] = 0.0;
			lm++;
		}
	} else {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = alm_im(shtns, im);
			l=im*MRES;
			mx[2*lm] = 0.0;
			while(++l <= LMAX) {
				a_1 = 1.0 / al[1];
				mx[2*lm+1] = a_1;		// specific to orthonormal.
				mx[2*lm+2] = a_1;
				lm++;	al+=2;
			}
			mx[2*lm +1] = 0.0;
			lm++;
		}
	}
}

/// fill mx with the coefficients of operator sin(theta).d/dtheta
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// stdq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void st_dt_matrix(shtns_cfg shtns, double* mx)
{
	mul_ct_matrix(shtns, mx);
	for (int lm=0; lm<NLM; lm++) {
		mx[2*lm]   *=   shtns->li[lm] - 1;
		mx[2*lm+1] *= -(shtns->li[lm] + 2);
	}
}

/// Multiplication of Qlm by a matrix involving l+1 and l-1 only.
/// The result is stored in Rlm, which MUST be different from Qlm.
/// mx is an array of 2*NLM values as returned by \ref mul_ct_matrix or \ref st_dt_matrix
/// compute: Rlm[lm] = mx[2*lm] * Qlm[lm-1] + mx[2*lm+1] * Qlm[lm+1];
void SH_mul_mx(shtns_cfg shtns, double* mx, cplx *Qlm, cplx *Rlm)
{
	long int nlmlim, lm;
	v2d* vq = (v2d*) Qlm;
	v2d* vr = (v2d*) Rlm;
	nlmlim = NLM-1;
	lm = 0;
		s2d mxu = vdup(mx[1]);
		vr[0] = mxu*vq[1];
	for (lm=1; lm<nlmlim; lm++) {
		s2d mxl = vdup(mx[2*lm]);		s2d mxu = vdup(mx[2*lm+1]);
		vr[lm] = mxl*vq[lm-1] + mxu*vq[lm+1];
	}
	lm = nlmlim;
		s2d mxl = vdup(mx[2*lm]);
		vr[lm] = mxl*vq[lm-1];
}

//@}

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX

/** \addtogroup local Local and partial evaluation of SH fields.
 * These do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
 * These functions are not optimized and can be relatively slow, but they provide good
 * reference implemenation for the transforms.
*/
//@{

/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
double SH_to_point(shtns_cfg shtns, cplx *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	double vr0, vr1;
	long int l,m,im;

	vr0 = 0.0;		vr1 = 0.0;
	m=0;	im=0;
		legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
		for (l=m; l<LTR; l+=2) {
			vr0 += yl[l] * creal( Qlm[l] );
			vr1 += yl[l+1] * creal( Qlm[l+1] );
		}
		if (l==LTR) {
			vr0 += yl[l] * creal( Qlm[l] );
		}
		vr0 += vr1;
	if (MTR>0) {
		cplx eip, eimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		im = 1;  do {
			m = im*MRES;
			legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
			v2d* Ql = (v2d*) &Qlm[LiM(shtns, 0,im)];	// virtual pointer for l=0 and im
			v2d vrm0 = vdup(0.0);		v2d vrm1 = vdup(0.0);
			for (l=m; l<LTR; l+=2) {
				vrm0 += vdup(yl[l]) * Ql[l];
				vrm1 += vdup(yl[l+1]) * Ql[l+1];
			}
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;			// not so accurate, but it should be enough for rendering uses.
			vrm0 += vrm1;
			if (l==LTR) {
				vrm0 += vdup(yl[l]) * Ql[l];
			}
			vr0 += vcplx_real(vrm0)*creal(eimp) - vcplx_imag(vrm0)*cimag(eimp);
		} while(++im <= MTR);
	}
	return vr0;
}

void SH_to_grad_point(shtns_cfg shtns, cplx *DrSlm, cplx *Slm, double cost, double phi,
					   double *gr, double *gt, double *gp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( DrSlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
		}
	if (MTR>0) {
		cplx eip, eimp, imeimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;		// not so accurate, but it should be enough for rendering uses.
			imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &DrSlm[l];		v2d* Sl = (v2d*) &Slm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dsdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dsdp += vdup(yl[l]) * Sl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);			// dS/dr
			vtt += vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp);		// dS/dt
			vpp += vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp);	// + I.m/sint *S
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*gr = vr0;	// Gr = dS/dr
	*gt = vtt;	// Gt = dS/dt
	*gp = vpp;	// Gp = I.m/sint *S
}

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
void SHqst_to_point(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( Qlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
			vpp -= dtyl[l] * creal( Tlm[l] );
		}
	if (MTR>0) {
		cplx eip, eimp, imeimp;
		eip = cos(phi*MRES) + I*sin(phi*MRES);	eimp = 2.0;
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
//			eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			eimp *= eip;		// not so accurate, but it should be enough for rendering uses.
			imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &Qlm[l];	v2d* Sl = (v2d*) &Slm[l];	v2d* Tl = (v2d*) &Tlm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dtdt = vdup(0.0);
			v2d dsdp = vdup(0.0);		v2d dtdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dtdt += vdup(dtyl[l]) * Tl[l];
				dsdp += vdup(yl[l]) * Sl[l];
				dtdp += vdup(yl[l]) * Tl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);
			vtt += (vcplx_real(dtdp)*creal(imeimp) - vcplx_imag(dtdp)*cimag(imeimp))	// + I.m/sint *T
					+ (vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp));	// + dS/dt
			vpp += (vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp))	// + I.m/sint *S
					- (vcplx_real(dtdt)*creal(eimp) - vcplx_imag(dtdt)*cimag(eimp));	// - dT/dt
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*vr = vr0;
	*vt = vtt;	// Bt = I.m/sint *T  + dS/dt
	*vp = vpp;	// Bp = I.m/sint *S  - dT/dt
}
//@}
	
#undef LTR
#undef MTR



/*
	SYNTHESIS AT A GIVEN LATITUDE
	(does not require a previous call to shtns_set_grid)
*/

fftw_plan ifft_lat = NULL;		///< fftw plan for SHqst_to_lat
int nphi_lat = 0;			///< nphi of previous SHqst_to_lat
double* ylm_lat = NULL;
double* dylm_lat;
double ct_lat = 2.0;
double st_lat;

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr, vt, and vp arrays must have nphi+2 doubles allocated (fftw requirement).
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe.
/// \ingroup local
void SHqst_to_lat(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr)
{
	cplx vst, vtt, vsp, vtp, vrr;
	cplx *vrc, *vtc, *vpc;
	long int m, l, j;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	vrc = (cplx *) vr;
	vtc = (cplx *) vt;
	vpc = (cplx *) vp;

	if ((nphi != nphi_lat)||(ifft_lat == NULL)) {
		if (ifft_lat != NULL) fftw_destroy_plan(ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		nphi_lat = nphi;
	}
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* NLM*2);
		dylm_lat = ylm_lat + NLM;
	}
	if (cost != ct_lat) {		// don't recompute if same latitude (ie equatorial disc rendering)
		st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
		for (m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	for (m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;	vtc[m] = 0.0;	vpc[m] = 0.0;
	}
	j=0;
	m=0;
		vrr=0;	vtt=0;	vst=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
			vst += dylm_lat[j] * creal(Slm[j]);
			vtt += dylm_lat[j] * creal(Tlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
		vtc[m] =  vst;	// Vt =   dS/dt
		vpc[m] = -vtt;	// Vp = - dT/dt
	for (m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;	vtt=0;	vst=0;	vsp=0;	vtp=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
			vst += dylm_lat[j] * Slm[j];
			vtt += dylm_lat[j] * Tlm[j];
			vsp += ylm_lat[j] * Slm[j];
			vtp += ylm_lat[j] * Tlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
		vtc[m] = I*m*vtp + vst;	// Vt = I.m/sint *T  + dS/dt
		vpc[m] = I*m*vsp - vtt;	// Vp = I.m/sint *S  - dT/dt
	}
	fftw_execute_dft_c2r(ifft_lat,vrc,vr);
	fftw_execute_dft_c2r(ifft_lat,vtc,vt);
	fftw_execute_dft_c2r(ifft_lat,vpc,vp);
//	free(ylm_lat);
}

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr arrays must have nphi+2 doubles allocated (fftw requirement).
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe.
/// \ingroup local
void SH_to_lat(shtns_cfg shtns, cplx *Qlm, double cost,
					double *vr, int nphi, int ltr, int mtr)
{
	cplx vrr;
	cplx *vrc;
	long int m, l, j;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	vrc = (cplx *) vr;

	if ((nphi != nphi_lat)||(ifft_lat == NULL)) {
		if (ifft_lat != NULL) fftw_destroy_plan(ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		nphi_lat = nphi;
	}
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* NLM*2);
		dylm_lat = ylm_lat + NLM;
	}
	if (cost != ct_lat) {		// don't recompute if same latitude (ie equatorial disc rendering)
		st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
		for (m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	for (m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;
	}
	j=0;
	m=0;
		vrr=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
	for (m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;
		for(l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
	}
	fftw_execute_dft_c2r(ifft_lat,vrc,vr);
//	free(ylm_lat);
}


/// complex scalar transform.
/// in: complex spatial field.
/// out: alm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of (LMAX+1)^2 coefficients.
void spat_cplx_to_SH(shtns_cfg shtns, cplx *z, cplx *alm)
{
	long int nspat = shtns->nspat;
	double *re, *im;
	cplx *rlm, *ilm;

	if (MMAX != LMAX) shtns_runerr("complex SH requires lmax=mmax and mres=1.");

	// alloc temporary fields
	re = (double*) VMALLOC( 2*(nspat + NLM*2)*sizeof(double) );
	im = re + nspat;
	rlm = (cplx*) (re + 2*nspat);
	ilm = rlm + NLM;

	// split z into real and imag parts.
	for (int k=0; k<nspat; k++) {
		re[k] = creal(z[k]);		im[k] = cimag(z[k]);
	}

	// perform two real transforms:
	spat_to_SH(shtns, re, rlm);
	spat_to_SH(shtns, im, ilm);

	// combine into complex coefficients
	int ll = 0;
	int lm = 0;
	for (int l=0; l<=LMAX; l++) {
		ll += 2*l;		// ll = l*(l+1)
		alm[ll] = creal(rlm[lm]) + I*creal(ilm[lm]);		// m=0
		lm++;
	}
	for (int m=1; m<=MMAX; m++) {
		ll = (m-1)*m;
		for (int l=m; l<=LMAX; l++) {
			ll += 2*l;		// ll = l*(l+1)
			cplx rr = rlm[lm];
			cplx ii = ilm[lm];
			alm[ll+m] = rr + I*ii;			// m>0
			rr = conj(rr) + I*conj(ii);		// m<0, m even
			if (m&1) rr = -rr;				// m<0, m odd
			alm[ll-m] = rr;
			lm++;
		}
	}

	VFREE(re);
}

/// complex scalar transform.
/// in: alm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of (LMAX+1)^2 coefficients.
/// out: complex spatial field.
void SH_to_spat_cplx(shtns_cfg shtns, cplx *alm, cplx *z)
{
	long int nspat = shtns->nspat;
	double *re, *im;
	cplx *rlm, *ilm;

	if (MMAX != LMAX) shtns_runerr("complex SH requires lmax=mmax and mres=1.");

	// alloc temporary fields
	re = (double*) VMALLOC( 2*(nspat + NLM*2)*sizeof(double) );
	im = re + nspat;
	rlm = (cplx*) (re + 2*nspat);
	ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	int ll = 0;
	int lm = 0;
	for (int l=0; l<=LMAX; l++) {
		ll += 2*l;		// ll = l*(l+1)
		rlm[lm] = creal(alm[ll]);		// m=0
		ilm[lm] = cimag(alm[ll]);
		lm++;
	}
	double half_parity = 0.5;
	for (int m=1; m<=MMAX; m++) {
		ll = (m-1)*m;
		half_parity = -half_parity;		// (-1)^m * 0.5
		for (int l=m; l<=LMAX; l++) {
			ll += 2*l;		// ll = l*(l+1)
			cplx b = alm[ll-m] * half_parity;		// (-1)^m for m negative.
			cplx a = alm[ll+m] * 0.5;
			rlm[lm] = (conj(b) + a);		// real part
			ilm[lm] = (conj(b) - a)*I;		// imag part
			lm++;
		}
	}

	// perform two real transforms:
	SH_to_spat(shtns, rlm, re);
	SH_to_spat(shtns, ilm, im);

	// combine into z
	for (int k=0; k<nspat; k++)
		z[k] = re[k] + I*im[k];

	VFREE(re);
}

/*
void SH_to_spat_grad(shtns_cfg shtns, cplx *alm, double *gt, double *gp)
{
	double *mx;
	cplx *blm, *clm;
	
	blm = (cplx*) VMALLOC( 3*NLM*sizeof(cplx) );
	clm = blm + NLM;
	mx = (double*)(clm + NLM);

	st_dt_matrix(shtns, mx);
	SH_mul_mx(shtns, mx, alm, blm);
	int lm=0;
	for (int im=0; im<=MMAX; im++) {
		int m = im*MRES;
		for (int l=m; l<=LMAX; l++) {
			clm[lm] = alm[lm] * I*m;
			lm++;
		}
	}
	SH_to_spat(shtns, blm, gt);
	SH_to_spat(shtns, clm, gp);
	for (int ip=0; ip<NPHI; ip++) {
		for (int it=0; it<NLAT; it++) {
			gt[ip*NLAT+it] /= shtns->st[it];
			gp[ip*NLAT+it] /= shtns->st[it];
		}
	}
	VFREE(blm);
}
*/
