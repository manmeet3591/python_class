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
	#include "SHT/spat_to_SHst_mic.c"
	#include "SHT/SHst_to_spat_mic.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SH_mic.c"
	#include "SHT/SH_to_spat_mic.c"
	#include "SHT/spat_to_SHst_mic.c"
	#include "SHT/SHst_to_spat_mic.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SH_mic.c"
	#include "SHT/SH_to_spat_mic.c"
	#include "SHT/spat_to_SHst_mic.c"
	#include "SHT/SHst_to_spat_mic.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SH_mic.c"
	#include "SHT/SH_to_spat_mic.c"
	#undef NWAY
	#define NWAY 6
	#include "SHT/spat_to_SH_mic.c"
	#include "SHT/SH_to_spat_mic.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SH_mic.c"
	#include "SHT/SH_to_spat_mic.c"
	#undef NWAY

#define SHT_GRAD
	#define NWAY 1
	#include "SHT/SHs_to_spat_mic.c"
	#include "SHT/SHt_to_spat_mic.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/SHs_to_spat_mic.c"
	#include "SHT/SHt_to_spat_mic.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_mic.c"
	#include "SHT/SHt_to_spat_mic.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHs_to_spat_mic.c"
	#include "SHT/SHt_to_spat_mic.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
	#define NWAY 1
	#include "SHT/spat_to_SHqst_mic.c"
	#include "SHT/SHqst_to_spat_mic.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SHqst_mic.c"
	#include "SHT/SHqst_to_spat_mic.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SHqst_mic.c"
	#include "SHT/SHqst_to_spat_mic.c"
	#undef NWAY
#undef SHT_3COMP


void* ffly[6][SHT_NTYP] = {
	{ NULL, NULL, SHsphtor_to_spat_mic1_l, spat_to_SHsphtor_mic1_l,
		SHsph_to_spat_mic1_l, SHtor_to_spat_mic1_l, SHqst_to_spat_mic1_l, spat_to_SHqst_mic1_l },
	{ SH_to_spat_mic2_l, spat_to_SH_mic2_l, SHsphtor_to_spat_mic2_l, spat_to_SHsphtor_mic2_l,
		SHsph_to_spat_mic2_l, SHtor_to_spat_mic2_l, SHqst_to_spat_mic2_l, spat_to_SHqst_mic2_l },
	{ SH_to_spat_mic3_l, spat_to_SH_mic3_l, SHsphtor_to_spat_mic3_l, spat_to_SHsphtor_mic3_l,
		SHsph_to_spat_mic3_l, SHtor_to_spat_mic3_l, SHqst_to_spat_mic3_l, spat_to_SHqst_mic3_l },
	{ SH_to_spat_mic4_l, spat_to_SH_mic4_l, NULL, NULL,
		SHsph_to_spat_mic4_l, SHtor_to_spat_mic4_l, NULL, NULL },
	{ SH_to_spat_mic6_l, spat_to_SH_mic6_l, NULL, NULL,
		NULL, NULL, NULL, NULL },
	{ SH_to_spat_mic8_l, spat_to_SH_mic8_l, NULL, NULL,
		NULL, NULL, NULL, NULL }
};

void* ffly_m0[6][SHT_NTYP] = {
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
};

void* fomp[6][SHT_NTYP] = {
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
};

void* ffly_m[6][SHT_NTYP] = {
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
	{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
};
