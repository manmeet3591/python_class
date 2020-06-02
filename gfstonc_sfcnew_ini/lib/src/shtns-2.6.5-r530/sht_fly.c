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

#include "sht_private.h"

#define MTR MMAX
#define SHT_VAR_LTR

#define GEN(name,sfx) GLUE2(name,sfx)
#define GEN3(name,nw,sfx) GLUE3(name,nw,sfx)

// genaral case
#undef SUFFIX
#define SUFFIX _l

	#define NWAY 1
	#include "SHT/spat_to_SHst_fly.c"
	#include "SHT/SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#include "SHT/spat_to_SHst_fly.c"
	#include "SHT/SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#include "SHT/spat_to_SHst_fly.c"
	#include "SHT/SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#undef NWAY

#define SHT_GRAD
	#define NWAY 1
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
	#define NWAY 1
	#include "SHT/spat_to_SHqst_fly.c"
	#include "SHT/SHqst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SHqst_fly.c"
	#include "SHT/SHqst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SHqst_fly.c"
	#include "SHT/SHqst_to_spat_fly.c"
	#undef NWAY
#undef SHT_3COMP

// axisymmetric
#define SHT_AXISYM
#undef SUFFIX
#define SUFFIX _m0l

	#define NWAY 1
	#include "SHT/spat_to_SHst_fly.c"
	#include "SHT/SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#include "SHT/spat_to_SHst_fly.c"
	#include "SHT/SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#include "SHT/spat_to_SHst_fly.c"
	#include "SHT/SHst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SH_fly.c"
	#include "SHT/SH_to_spat_fly.c"
	#undef NWAY

#define SHT_GRAD
	#define NWAY 1
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHs_to_spat_fly.c"
	#include "SHT/SHt_to_spat_fly.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
	#define NWAY 1
	#include "SHT/spat_to_SHqst_fly.c"
	#include "SHT/SHqst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SHqst_fly.c"
	#include "SHT/SHqst_to_spat_fly.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SHqst_fly.c"
	#include "SHT/SHqst_to_spat_fly.c"
	#undef NWAY
#undef SHT_3COMP



void* ffly[6][SHT_NTYP] = {
	{ NULL, NULL, SHsphtor_to_spat_fly1_l, spat_to_SHsphtor_fly1_l,
		SHsph_to_spat_fly1_l, SHtor_to_spat_fly1_l, SHqst_to_spat_fly1_l, spat_to_SHqst_fly1_l },
	{ SH_to_spat_fly2_l, spat_to_SH_fly2_l, SHsphtor_to_spat_fly2_l, spat_to_SHsphtor_fly2_l,
		SHsph_to_spat_fly2_l, SHtor_to_spat_fly2_l, SHqst_to_spat_fly2_l, spat_to_SHqst_fly2_l },
	{ SH_to_spat_fly3_l, spat_to_SH_fly3_l, SHsphtor_to_spat_fly3_l, spat_to_SHsphtor_fly3_l,
		SHsph_to_spat_fly3_l, SHtor_to_spat_fly3_l, SHqst_to_spat_fly3_l, spat_to_SHqst_fly3_l },
	{ SH_to_spat_fly4_l, spat_to_SH_fly4_l, NULL, NULL,
		SHsph_to_spat_fly4_l, SHtor_to_spat_fly4_l, NULL, NULL },
	{ SH_to_spat_fly6_l, spat_to_SH_fly6_l, NULL, NULL,
		NULL, NULL, NULL, NULL },
	{ SH_to_spat_fly8_l, spat_to_SH_fly8_l, NULL, NULL,
		NULL, NULL, NULL, NULL }
};

void* ffly_m0[6][SHT_NTYP] = {
	{ NULL, NULL, SHsphtor_to_spat_fly1_m0l, spat_to_SHsphtor_fly1_m0l,
		SHsph_to_spat_fly1_m0l, SHtor_to_spat_fly1_m0l, SHqst_to_spat_fly1_m0l, spat_to_SHqst_fly1_m0l },
	{ SH_to_spat_fly2_m0l, spat_to_SH_fly2_m0l, SHsphtor_to_spat_fly2_m0l, spat_to_SHsphtor_fly2_m0l,
		SHsph_to_spat_fly2_m0l, SHtor_to_spat_fly2_m0l, SHqst_to_spat_fly2_m0l, spat_to_SHqst_fly2_m0l },
	{ SH_to_spat_fly3_m0l, spat_to_SH_fly3_m0l, SHsphtor_to_spat_fly3_m0l, spat_to_SHsphtor_fly3_m0l,
		SHsph_to_spat_fly3_m0l, SHtor_to_spat_fly3_m0l, SHqst_to_spat_fly3_m0l, spat_to_SHqst_fly3_m0l },
	{ SH_to_spat_fly4_m0l, spat_to_SH_fly4_m0l, NULL, NULL,
		SHsph_to_spat_fly4_m0l, SHtor_to_spat_fly4_m0l, NULL, NULL },
	{ SH_to_spat_fly6_m0l, spat_to_SH_fly6_m0l, NULL, NULL,
		NULL, NULL, NULL, NULL },
	{ SH_to_spat_fly8_m0l, spat_to_SH_fly8_m0l, NULL, NULL,
		NULL, NULL, NULL, NULL }
};

void* ffly_m[6][SHT_NTYP] = {
	{ NULL, NULL, SHsphtor_m_to_spat_fly1_l, spat_to_SHsphtor_m_fly1_l,
		SHsph_m_to_spat_fly1_l, SHtor_m_to_spat_fly1_l, SHqst_m_to_spat_fly1_l, spat_to_SHqst_m_fly1_l },
	{ SH_m_to_spat_fly2_l, spat_to_SH_m_fly2_l, SHsphtor_m_to_spat_fly2_l, spat_to_SHsphtor_m_fly2_l,
		SHsph_m_to_spat_fly2_l, SHtor_m_to_spat_fly2_l, SHqst_m_to_spat_fly2_l, spat_to_SHqst_m_fly2_l },
	{ SH_m_to_spat_fly3_l, spat_to_SH_m_fly3_l, SHsphtor_m_to_spat_fly3_l, spat_to_SHsphtor_m_fly3_l,
		SHsph_m_to_spat_fly3_l, SHtor_m_to_spat_fly3_l, SHqst_m_to_spat_fly3_l, spat_to_SHqst_m_fly3_l },
	{ SH_m_to_spat_fly4_l, spat_to_SH_m_fly4_l, NULL, NULL,
		SHsph_m_to_spat_fly4_l, SHtor_m_to_spat_fly4_l, NULL, NULL },
	{ SH_m_to_spat_fly6_l, spat_to_SH_m_fly6_l, NULL, NULL,
		NULL, NULL, NULL, NULL },
	{ SH_m_to_spat_fly8_l, spat_to_SH_m_fly8_l, NULL, NULL,
		NULL, NULL, NULL, NULL }
};
