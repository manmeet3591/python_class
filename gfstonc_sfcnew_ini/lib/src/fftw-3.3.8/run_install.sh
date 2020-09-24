#!/bin/bash
./configure --prefix=/home/manmeet/Documents/python_class/gfstonc_sfcnew_ini/lib --enable-openmp --enable-shared --enable-avx i
make
make install
export FFTW3_LIBDIR=/home/manmeet/Documents/python_class/gfstonc_sfcnew_ini/lib/lib
