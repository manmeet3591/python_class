**SHTns is a high performance library for Spherical Harmonic Transform written in C,
aimed at numerical simulation (fluid flows, mhd, ...) in spherical geometries.**

Copyright (c) 2010-2014 Centre National de la Recherche Scientifique.
written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
SHTns is distributed under the open source [CeCILL License](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html)
(GPL compatible) located in the LICENSE file.

FEATURES:
---------

- **blazingly fast**.
- both **scalar and vector transforms**.
- backward and forward (synthesis and analysis) functions.
- flexible truncation (degree, order, azimuthal periodicity).
- spatial data can be stored in latitude-major or longitude-major arrays.
- various conventions (normalization and Condon-Shortley phase).
- can be used from **Fortran, c/c++, and Python** programs.
- a highly efficient Gauss algorithm working with Gauss nodes (based on
  Gauss-Legendre quadrature).
- support for SSE2, SSE3 and **AVX** vectorization, as well as Xeon Phi and
  Blue Gene/Q.
- **parallel transforms with OpenMP** (for Gauss grid only).
- an algorithm using DCT for regular nodes (generalized Fejer quadrature).
- synthesis (inverse transform) at any coordinate (not constrained to a grid).
- ability to choose the optimal spatial sizes for a given spectral truncation.
- **on-the-fly transforms** : saving memory and bandwidth, they are even faster
  on modern architectures.
- accurate up to spherical harmonic degree l=16383 (at least).
- rotation functions to rotate spherical harmonics (beta).
- special spectral operator functions that do not require a transform
  (multiply by cos(theta)...).
- scalar transforms for complex spatial fields.
- SHT at fixed m (without fft, aka Legendre transform - beta).


INSTALL:
--------

Briefly, the shell commands `./configure; make; make install` should
configure, build, and install this package. `./configure --help` will
list available options (among which `--enable-openmp` and `--enable-python`).
However, in order to get the best performance, it is highly recommended to
compile and install the FFTW library yourself, because many distributions
include a non-optimized FFTW library.

DOCUMENTATION:
--------------

- On-line doc is available: <http://users.isterre.fr/nschaeff/SHTns/>
- You can build it locally: Run `make docs` to generate documentation
  (requires doxygen). 
  Then browse the html documentation starting with `doc/html/index.html`
- A related research paper has been published:
  [Efficient Spherical Harmonic Transforms aimed at pseudo-spectral numerical simulations](http://dx.doi.org/10.1002/ggge.20071),
  also [available from arXiv](http://arxiv.org/abs/1202.6522).
- If you use SHTns for research work, please **cite this paper**:

        @article {shtns,
          author = {Schaeffer, Nathanael},
          title = {Efficient spherical harmonic transforms aimed at
          pseudospectral numerical simulations},
          journal = {Geochemistry, Geophysics, Geosystems},
          doi = {10.1002/ggge.20071},
          volume = {14}, number = {3}, pages = {751--758},
          year = {2013},
        }

CHANGE LOG:
-----------

* v2.6.5  (24 Aug 2015)
	- critical bugfix (multiple transforms of different lmax failed sometimes).
	- faster DCT initialization when using multiple threads.

* v2.6.4  (25 Jul 2015)
	- a critical bugfix (segfault of fixed-m tranforms when mmax=0).
	- fixed-m Legendre transforms added to Fortran API.

* v2.6.3  (9 Mar 2015)
	- better default compilation flags for icc.
	- complex transforms added to Fortran API (thanks to Bertrand Putigny)

* v2.6.2  (30 Dec 2014)
	- fix regression: Schmidt normalized analysis failed since v2.6 in some cases.

* v2.6.1  (17 Dec 2014)
	- new functions in python interface to control console output.
	- fix: `spat_cplx_to_SH()` and `SH_to_spat_cplx()` were missing
	  a (-1)^m for m<0 [issue #16].
	- fix: segfault in `spat_to_SH_ml()` [issue #15].
	- fix a few compilation issues.

* v2.6  (24 Oct 2014)
	- support for IBM Blue Gene/Q (QPX) with [bgclang](http://trac.alcf.anl.gov/projects/llvm-bgq).
          Configure with `./configure --enable-many-core CC=bgclang`
	- new beta feature: SHT at fixed m (aka Legendre transform).
	- faster initialization with OpenMP.
	- fix: in python, a rare coredump now correctly raises an exception.
	- fix: a few compilation problems.

* v2.5  (13 Mar 2014)
	- new experimental support for Intel Xeon Phi (MIC) in native mode,
	  with contributions from Vincent Boulos (Bull). For good performance
	  icc 14 is required. Configure with `./configure --enable-mic CC=icc`
	- fftw3.h included for easier compilation.
	- fix: obey `OMP_NUM_THREADS` environement variable
	- fix: failure of fly analysis with some special (rare) sizes.
	- add missing `shtns_print_cfg()` to Fortran interface
	- new save/restore plan feature for bit-level reproducibility

* v2.4.1  (18 Sep 2013)
	- performance improvement: analysis with `SHT_PHI_CONTIGUOUS` is now
	  on par with synthesis (or better), even for large transforms.

* v2.4  (5 Aug 2013)
	- new scalar transforms for complex spatial fields: `SH_to_spat_cplx()`
	  and `spat_cplx_to_SH()`.
	- new `shtns_verbose()` function to control output during initialization.
	- better MKL support (including multi-thread). Warning: MKL's FFTW
	  interface is not thread safe, SHTns can't be called from multiple
	  threads if compiled with MKL.
	- fix compatibility with c++ std::complex.
	- new shallow water simulation example in examples/

* v2.3.1  (10 Apr 2013)
	- OpenMP library is now installed as `libshtns_omp.a`.
	- fix detection of OpenMP mutlithreaded FFTW.
	- new configure option `--enable-mkl` to use the FFT of the MKL
	  library instead of FFTW (lower performance expected).
	- `time_SHT` can be compiled on MacOSX and uses less memory.
	- new `SH_to_lat()` function.
	- a few other minor improvements and fixes.

* v2.3  (3 Oct 2012)
	- added `mi` member in `shtns_info` structure (ABI change).
	- added function to access the Gauss nodes.
	- added support for special operators in spectral space (multiplication
	  by cos(theta) and sin(theta).d/dtheta for instance).
	- shtns.h is now compatible with C++.
	- better python interface for rotations.
	- performance improvement for OpenMP code without `fftw3_omp`.
	- slightly faster `SH_to_point()` [5%] and `SHqst_to_point()` [20%].
	- bugfix: in some rare cases, OpenMP code freed unallocated memory.
	- bugfix: fixed python interface compilation with clang.

* v2.2.4  (25 Jun 2012)
	- the previous critical bugfix had not been applied to parallel OpenMP
	  transforms.

* v2.2.3  (24 Jun 2012)
	- critical bugfix: `SHtor_to_spat()` and `SHsph_to_spat()` gave wrong results
	  for mmax>0 with on-the-fly transforms.
	- minor bugfix in Python interface.

* v2.2.2  (21 Jun 2012)
	- better Python interface: using `synth()` and `analys()` methods.
	- bugfix in build system: can now compile python extension without openmp.

* v2.2.1  (21 May 2012)
	- slightly faster parallel transforms.
	- better Python interface: decent error handling and keyword argument support.
	- changes to Python interface: grid defaults to `SHT_PHI_CONTIGUOUS`, 
	  `set_grid_auto()` removed.
	- bugfix: default compilation with FFTW 3.0 to avoid "bad Gauss points" error.
	- bugfix: correct alignement of gauss weights in 32 bit systems to avoid
	  segfaults.
	- new ./configure script for easier configuration and compilation.

* v2.2  (23 Apr 2012)
	- parallel transforms with OpenMP (for Gauss grid, significant benefit
	  for l>=127).

* v2.1  (8 Mar 2012)
	- support for huge spherical harmonic degree (tested up to l>43600).
	- speed improvements, especially for large transforms.
	- compilation with FFTW v3.0 or more is now possible through a 
	  configuration option (see `sht_config.h`)

* v2.0  (9 Feb 2012)
	- support for AVX instruction set (almost x2 speed-up on Sandy-Bridge
	  processors).
	- allow multiple transforms with different sizes, normalizations and
	  grids (C interface only).
	- changes to C interface : most functions now require a handle to identify
	  the transform. (Fortran interface unchanged)
	- transforms are accurate up to spherical harmonic degree l=2700 (at least).
	- lots of small improvements, speed-ups and a few bug fixes.
	- requires FFTW v3.3.
	- better Python interface using NumPy arrays (beta).
	- rotation functions to rotate spherical harmonics (beta).

* v1.5  (4 May 2011)
	- on-the-fly transforms which do not require huge matrices : save memory
	  and bandwidth, and can be faster on some architecture.
	- runtime selection of fastest algorithm, including on-the-fly transforms.
	- transforms are accurate up to spherical harmonic degree l=2045 (at least).
	- fix a bug that lead to wrong results for `SHtor_to_spat` and `SHsph_to_spat`.
	- a bunch of minor improvements, optimizations and fixes.

* v1.0  (9 June 2010)
	- initial release for C/C++ and Fortran under CeCILL licence (GPL compatible).
	- scalar and vector, forward and backward transforms.
	- support several normalization conventions.
	- transforms are accurate up to spherical harmonic degree l=1300 (at least).
	- flexible truncation and spatial sizes.
	- support spatial data stored in latitude-major or longitude-major arrays.
	- regular grid (with DCT acceleration) or Gauss grid (highly optimized).
	- SSE2 vectorization.
	- synthesis at any coordinate (not constrained to grid).
	- can choose the optimal spatial size for a given spectral truncation.
	- requires FFTW 3.0.
