#!/bin/bash
#source activate gfs27
./configure --prefix=/home/manmeet/Documents/python_class/gfstonc_sfcnew_ini/lib --enable-python
make
make install-lib
export SHTNS_LIBDIR=/home/manmeet/Documents/python_class/gfstonc_sfcnew_ini/lib/lib
cp libshtns_omp.a /home/manmeet/Documents/python_class/gfstonc_sfcnew_ini/lib/lib/libshtns.a
python setup.py install
