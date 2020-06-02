## path prefix for make install (installs to $(PREFIX)/lib and $(PREFIX)/include)
PREFIX=/home/manmeet/Documents/gfstonc/lib
SED=/bin/sed
CC=gcc -march=native
CFLAGS=-O2 -I/home/manmeet/Documents/gfstonc/lib/include -L/home/manmeet/Documents/gfstonc/lib/lib -ffast-math -fomit-frame-pointer -std=gnu99 -fopenmp -fPIC 
LIBS=-lfftw3_omp -lfftw3 -lrt -lm 
PYTHON=python


## global options for gcc
## there should be -ffast-math or at least -fcx-limited-range to produce fast code.
go= $(CFLAGS) -D_GNU_SOURCE


# compiler command with options for the sht codelets
shtcc = $(CC) $(go) -O3
# compiler command with options for other source (initialization, ...)
cc = $(CC) $(go) -O2

hfiles = sht_private.h sht_config.h shtns.h
objs = sht_init.o sht_fly.o sht_omp.o sht_mem.o sht_dct.o
libname = libshtns_omp.a

default : _shtns.so

install : install-py

libshtns_omp.a : Makefile $(objs)
	ar rcs libshtns_omp.a $(objs)
	@echo " "
	@cat COPYRIGHT
	@echo " "
	@echo "**** link with :  -lfftw3_omp -lfftw3 -lrt -lm "

install-lib : $(libname)
	@mkdir -p $(PREFIX)/lib/
	@mkdir -p $(PREFIX)/include/
	cp $(libname) $(PREFIX)/lib/$(libname)
	cp shtns.h $(PREFIX)/include/
	cp shtns.f $(PREFIX)/include/
	@echo " "
	@cat COPYRIGHT
	@echo " "
	@echo "**** link with : -L$(PREFIX)/lib  -lfftw3_omp -lfftw3 -lrt -lm "

# codelets :
SHT/SH_to_spat.c : SHT/hyb_SH_to_spat.gen.c
	$(MAKE) SH_to_spat.c -C SHT SED=$(SED)
SHT/spat_to_SH.c : SHT/hyb_spat_to_SH.gen.c
	$(MAKE) spat_to_SH.c -C SHT SED=$(SED)
SHT/SH_to_spat_fly.c : SHT/fly_SH_to_spat.gen.c
	$(MAKE) SH_to_spat_fly.c -C SHT SFX=fly SED=$(SED)
SHT/spat_to_SH_fly.c : SHT/fly_spat_to_SH.gen.c
	$(MAKE) spat_to_SH_fly.c -C SHT SFX=fly SED=$(SED)
SHT/SH_to_spat_omp.c : SHT/omp_SH_to_spat.gen.c
	$(MAKE) SH_to_spat_omp.c -C SHT SFX=omp SED=$(SED)
SHT/spat_to_SH_omp.c : SHT/omp_spat_to_SH.gen.c
	$(MAKE) spat_to_SH_omp.c -C SHT SFX=omp SED=$(SED)
SHT/SH_to_spat_mic.c : SHT/mic_SH_to_spat.gen.c
	$(MAKE) SH_to_spat_mic.c -C SHT SFX=mic SED=$(SED)
SHT/spat_to_SH_mic.c : SHT/mic_spat_to_SH.gen.c
	$(MAKE) spat_to_SH_mic.c -C SHT SFX=mic SED=$(SED)

# objects :
sht_init.o : sht_init.c Makefile sht_legendre.c sht_func.c sht_com.c $(hfiles)
	$(cc) -c $< -o $@

sht_mem.o : sht_mem.c Makefile $(hfiles) SHT/SH_to_spat.c SHT/spat_to_SH.c
	$(shtcc) -c $< -o $@
sht_dct.o : sht_dct.c Makefile $(hfiles) SHT/SH_to_spat.c SHT/spat_to_SH.c
	$(shtcc) -c $< -o $@
sht_fly.o : sht_fly.c Makefile $(hfiles) SHT/SH_to_spat_fly.c SHT/spat_to_SH_fly.c
	$(shtcc) -c $< -o $@
sht_omp.o : sht_omp.c Makefile $(hfiles) SHT/SH_to_spat_omp.c SHT/spat_to_SH_omp.c
	$(shtcc) -c $< -o $@
sht_mic.o : sht_mic.c Makefile $(hfiles) SHT/SH_to_spat_mic.c SHT/spat_to_SH_mic.c
	$(cc) -c $< -o $@


# programs & Examples
time_SHT : time_SHT.c $(libname) shtns.h
	$(cc) time_SHT.c -I. ./$(libname) $(LIBS) -o time_SHT

test_rot : examples/test_rot.c $(libname) shtns.h
	$(cc) examples/test_rot.c -I. ./$(libname) $(LIBS) -o test_rot

SHT_example : examples/SHT_example.c $(libname) shtns.h
	$(cc) examples/SHT_example.c -I. ./$(libname) $(LIBS) -o SHT_example

SHT_fort_ex : examples/SHT_example.f $(libname) shtns.f
	gfortran -L$(PREFIX)/lib -I$(PREFIX)/include -I. -fdefault-real-8 examples/SHT_example.f ./$(libname) $(LIBS) -lc -o SHT_fort_ex

#documentation :
docs :
	( cat doc/doxygen.conf ; echo "PROJECT_NUMBER=2.6.5" ) | doxygen -

clean :
	$(MAKE) clean -C SHT
	rm -f *.o *.a
	rm -rf doc/html
	rm -rf doc/latex
	rm -rf build


# build a python interface using SWIG.
# use it with "from shtns import *" in a python program/shell
_shtns.so : shtns_numpy_wrap.c Makefile $(hfiles) $(objs)
	$(PYTHON) setup.py build
	@echo "*****************************************************************************************"
	@echo "*** Now, for system wide install, run 'make install' as root.                           *"
	@echo "*** Otherwise run '$(PYTHON) setup.py install --user' to install for current user only. *"
	@echo "*****************************************************************************************"

# generate python and c glue code with SWIG.
shtns_numpy_wrap.c : shtns_numpy.i sht_private.h shtns.h
	swig -python shtns_numpy.i

install-py : _shtns.so
	$(PYTHON) setup.py install

.PHONY : install install-py install-lib clean docs


#fftw compiling options :
#-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -fno-schedule-insns -fno-web -fno-loop-optimize --param inline-unit-growth=1000 --param large-function-growth=1000

