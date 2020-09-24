#!/bin/bash
./configure --prefix=/home/manmeet/Documents/gfstonc/lib --enable-openmp --enable-shared --enable-avx 
make
make install
export FFTW3_LIBDIR=/home/manmeet/Documents/gfstonc/lib/lib
