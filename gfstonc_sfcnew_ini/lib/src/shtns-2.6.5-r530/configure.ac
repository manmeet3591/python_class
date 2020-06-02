#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.62])
AC_INIT([SHTns],[2.6.5],[],[shtns],[https://bitbucket.org/nschaeff/shtns])
AC_LANG([C])
AC_CONFIG_SRCDIR([sht_init.c])
AC_CONFIG_HEADERS([sht_config.h])
target="libshtns.a"		# by default, build the library.
objs="sht_fly.o"
install="install-lib"	# by default, install the library.

# optional variables :
AC_ARG_VAR(PYTHON, [the python interpreter (defaults to 'python')])
# optional features with --enable-XXX
AC_ARG_ENABLE([verbose],
	AS_HELP_STRING([--enable-verbose=0-3], [define verbosity level for the library: 0=silent, 1=default, 2=debug, 3=full]),
	[],[enable_verbose=1])
AC_ARG_ENABLE([openmp],
	AS_HELP_STRING([--enable-openmp], [Enable multi-threading with OpenMP]), [], [enable_openmp=default])
AC_ARG_ENABLE([mkl],
	AS_HELP_STRING([--enable-mkl], [Try to link with intel MKL instead of FFTW (can be slightly faster, but NOT THREAD SAFE)]), [], [enable_mkl=default])
AC_ARG_ENABLE([mic],
	AS_HELP_STRING([--enable-mic], [Enable compilation for Xeon Phi]), [], [enable_mic=no])
AC_ARG_ENABLE([many-core],
	AS_HELP_STRING([--enable-many-core], [Enable optimizations for many-core systems (Blue Gene/Q, Xeon Phi,...) and remove single-threaded transforms]), [], [enable_many_core=no])
AC_ARG_ENABLE([python],
	AS_HELP_STRING([--enable-python], [Build Python interface]), [], [enable_python=no])
AC_ARG_ENABLE([long-double],
	AS_HELP_STRING([--enable-long-double], [Use long double during initialization for (maybe) slightly better precision]))
AC_ARG_ENABLE([f77],
	AS_HELP_STRING([--disable-f77], [Do not include F77 wrapper to call SHTns library from Fortran]))
AC_ARG_ENABLE([simd],
	AS_HELP_STRING([--disable-simd], [Do not use vector extensions (SSE2, AVX or MIC)]))
AC_ARG_ENABLE([dct],
	AS_HELP_STRING([--disable-dct], [Do not use DCT, removes support for regular grids.]))
AC_ARG_ENABLE([mem],
	AS_HELP_STRING([--disable-mem], [Do not use matrix precomputed and stored in memory. Implies --disable-dct.]))

dnl Sanitize $prefix. Autoconf does this by itself, but so late in the
dnl generated configure script that the expansion does not occur until
dnl after our eval magic below.
AS_IF([test "$prefix" = "NONE"],[prefix=$ac_default_prefix])

# Checks for programs.
if test "x$CFLAGS" = "x"; then
	CFLAGS="-O2"
fi
CFLAGS="$CFLAGS -I$prefix/include -L$prefix/lib"
AC_PROG_CC
AC_PROG_SED
if test "$SED" = :; then
	AC_MSG_ERROR([sed program required.])
fi

# define macrao AX_CHECK_COMPILE_FLAG
AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_IF([test x"AS_VAR_GET(CACHEVAR)" = xyes],
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

AS_IF([test "x$enable_mic" != "xyes"], [
	# add gcc compile options if supported.
	AX_CHECK_COMPILE_FLAG([-march=native],[CC="$CC -march=native"], [
		AX_CHECK_COMPILE_FLAG([-mtune=native],[CC="$CC -mtune=native"])
	])
  ],[
	enable_many_core=yes		# optimize for many-core architecture
	# MIC native cross-compiling
	AX_CHECK_COMPILE_FLAG([-mmic],[CC="$CC -mmic"], [AC_MSG_ERROR(["Xeon Phi not supported."])])
	# With MIC, enable mkl by default.
	AS_IF([test "x$enable_mkl" = "xdefault"], [enable_mkl=yes])
])
AX_CHECK_COMPILE_FLAG([-ffast-math],[CFLAGS="$CFLAGS -ffast-math"])
AX_CHECK_COMPILE_FLAG([-fomit-frame-pointer],[CFLAGS="$CFLAGS -fomit-frame-pointer"])
AX_CHECK_COMPILE_FLAG([-std=gnu99],[CFLAGS="$CFLAGS -std=gnu99"])

# with intel compiler, we should not use -O3 (bad performance with icc14 and icc15)
AC_EGREP_CPP(icc,
  [#ifdef __INTEL_COMPILER
   icc
   #endif
  ], [shtcc_flags="-O2"], [shtcc_flags="-O3"])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h stdio.h string.h math.h complex.h])

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,AC_MSG_ERROR([math library not found.]))
# this library is used by time_SHT only (and not required at all on MacOSX)
AC_CHECK_LIB([rt],[clock_gettime],,)

# With Python, enable openmp by default.
AS_IF([test "x$enable_python" != "xno"], [
	AS_IF([test "x$enable_openmp" = "xdefault"], [enable_openmp=yes])
])

# On many-core systems, enable openmp by default, and disable precomputed matrix support.
AS_IF([test "x$enable_many_core" != "xno"], [
  	enable_mem=no
	AS_IF([test "x$enable_openmp" = "xdefault"], [enable_openmp=yes])
])

AS_IF([test "x$enable_openmp" = "xyes"], [
	AC_OPENMP
	CFLAGS="$CFLAGS $OPENMP_CFLAGS"
	dnl AC_CHECK_HEADERS([omp.h])
	AS_IF([test "x$enable_mic" != "xno"], [
		target="libshtns_mic.a"
		objs="sht_mic.o"
	],[
		target="libshtns_omp.a"
		objs="$objs sht_omp.o"
		AS_IF([test "x$enable_many_core" != "xno"], [objs="sht_mic.o"])
	])
])
libname=$target		# set the libname befor checking for python.


# Check for FFTW (FFTW preferred, MKL optional)
AS_IF([test "x$enable_mkl" = "xyes"], [
  # optionally use MKL instead of fftw3
  AC_DEFINE(HAVE_LIBFFTW3_OMP, 1, [MKL has the multi-thread fftw interface.])
  enable_dct=no		# MKL does not support dct (june 2013).
  AX_CHECK_COMPILE_FLAG([-mkl], [
	LIBS="-mkl $LIBS"	# easy way !
	AC_CHECK_FUNC(fftw_plan_many_dft, [AC_MSG_NOTICE(["FFTW interface found in MKL, link with -mkl"])],
		[AC_MSG_ERROR(["FFTW interface not found in MKL libraries (-mkl)."])])
  ],[
	AS_IF([test "x$enable_openmp" = "xyes"], [
		# multi-threaded MKL libs
		#mkl_libs="-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread"
		mkl_libs="-Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group"
	],[
		# sequential MKL libs
		mkl_libs="-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group"
	])
	AC_CHECK_LIB(mkl_intel_lp64, fftw_plan_many_dft, [LIBS="$mkl_libs $LIBS"],
		AC_MSG_ERROR(["FFTW interface not found in MKL libraries ($mkl_libs)."]),
	["$mkl_libs"])
	AC_MSG_NOTICE(["FFTW interface found in MKL, link with $mkl_libs"])
  ])
],[
	# fftw3
	AC_CHECK_LIB([fftw3],[fftw_plan_many_dft],,AC_MSG_ERROR([FFTW3 library required.]))
	AS_IF([test "x$enable_openmp" = "xyes"], [
		AC_CHECK_LIB([fftw3_omp], [fftw_init_threads],,AC_MSG_WARN([FFTW3 does not support OpenMP. Bad multi-threaded performance expected.]))
	])
])


# Checks related to Python and NumPy paths:
AS_IF([test "x$enable_python" != "xno"], [
	AS_IF([test "x$PYTHON" = "x"], [
		AC_MSG_CHECKING(for python with numpy package)
		py_test="python python2 python3"
	],[
		AC_MSG_CHECKING(if $PYTHON has numpy package)
		py_test="$PYTHON"
	])
	for py in $py_test
	do :
		numpy_inc=`$py -c "from numpy import get_include; print(get_include())" 2>/dev/null`
		AS_IF([test "x$numpy_inc" != "x"],[break])
	done
	AS_IF([test "x$numpy_inc" = "x"], [
		AC_MSG_RESULT(no)
		AC_MSG_ERROR([NumPy package is required for the python interface.])
	],[
		AC_MSG_RESULT($py)
		PYTHON=$py
	])
	python_inc=`$PYTHON -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())"  2>/dev/null`
	echo "  python include path = $python_inc"
	echo "  numpy include path = $numpy_inc"
	AC_SUBST([numpy_inc])
	AC_SUBST([python_inc])
	CFLAGS="$CFLAGS -fpic"	# required compile flag for python extensions.
	target="_shtns.so"		# build python extension instead of C library.
	install="install-py"	# install python extension instead of C library.
])

# Checks related to long double
AS_IF([test "x$enable_long_double" = "xyes"], [
	AC_TYPE_LONG_DOUBLE_WIDER])

# Disable Fortran interface ?
AS_IF([test "x$enable_f77" != "xno"], [
	AC_DEFINE([SHT_F77_API],[1],[Compile the Fortran API])
])

# Disable SIMD ?
AS_IF([test "x$enable_simd" != "xno"], [
	AC_DEFINE([_GCC_VEC_],[1],[I compile with GCC 4 or ICC 14 or later, and I would like fast vectorized code (if SSE2, AVX or MIC is supported) !])
])

# Disable MEM ?
AS_IF([test "x$enable_mem" != "xno"], [
	AC_DEFINE([SHTNS_MEM],[1],[Include algorithms using precomputed matrix stored in memory)])
	objs="$objs sht_mem.o"		# compile dct transforms
	# Disable DCT ?
	AS_IF([test "x$enable_dct" != "xno"], [
		AC_DEFINE([SHTNS_DCT],[1],[Include DCT support (Fejer quadrature on regular grids)])
		objs="$objs sht_dct.o"		# compile dct transforms
	])
])

# Verbosity setting
AS_IF([test "x$enable_verbose" == "xno"], [enable_verbose=0],
	[test "x$enable_verbose" == "xyes"], [enable_verbose=1])
AC_DEFINE_UNQUOTED([SHT_VERBOSE],$enable_verbose,[0:no output, 1:output info to stdout, 2:more output (debug info), 3:also print fftw plans.])
echo "  verbose level = $enable_verbose"

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([fftw_cost])

### for cycle.h ###
   AC_C_INLINE
   AC_HEADER_TIME
   AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])

   AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])

   AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])

   dnl Cray UNICOS _rtc() (real-time clock) intrinsic
   AC_MSG_CHECKING([for _rtc intrinsic])
   rtc_ok=yes
   AC_TRY_LINK([#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
   AC_MSG_RESULT($rtc_ok)
### end cycle.h ###

echo "prefix=$prefix"
AC_SUBST([target])
AC_SUBST([install])
AC_SUBST([objs])
AC_SUBST([libname])
AC_SUBST([shtcc_flags])
AC_CONFIG_FILES([Makefile setup.py])
AC_OUTPUT
