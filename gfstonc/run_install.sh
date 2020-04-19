#!/bin/bash
cd src

gfortran -O3 -march=native -O3 -fPIC -c kinds.f90
gfortran -O3 -march=native -O3 -fPIC -c sfcio_module.f90
gfortran -O3 -march=native -O3 -fPIC -c shtns.f90
#gfortran -O3 -march=native -fPIC -c nemsio_module.f90
#gfortran -O3 -march=native -fPIC -c nemsio_openclose.f90
#gfortran -O3 -march=native -fPIC -c nemsio_read.f90


python -m numpy.f2py read_sfc.f90 -m _read_sfc -h _read_sfc.pyf --overwrite-signature
python -m numpy.f2py -c _read_sfc.pyf read_sfc.f90

python -m numpy.f2py read_sigma_spec.f90 -m _read_sigma_spec -h _read_sigma_spec.pyf --overwrite-signature
python -m numpy.f2py -c _read_sigma_spec.pyf read_sigma_spec.f90

python -m numpy.f2py write_sfc.f90 -m _write_sfc -h _write_sfc.pyf --overwrite-signature
python -m numpy.f2py -c _write_sfc.pyf write_sfc.f90

#python -m numpy.f2py read_sfcflx_nemsio.f90 -m _read_sfcflx_nemsio -h _read_sfcflx_nemsio.pyf --overwrite-signature
#python -m numpy.f2py -c _read_sfcflx_nemsio.pyf read_sfcflx_nemsio.f90
#
#python -m numpy.f2py read_sigma_nemsio.f90 -m _read_sigma_nemsio -h _read_sigma_nemsio.pyf --overwrite-signature
#python -m numpy.f2py -c _read_sigma_nemsio.pyf read_sigma_nemsio.f90

cd -
rm -rf build/*
python setup.py build
python setup.py install
