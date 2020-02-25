from numpy.distutils.core  import setup, Extension
import os, sys, subprocess

# build fortran libraries if they do not yet exist.
if not os.path.isfile('src/libw3nco_d.a'):
    strg = 'cd src/w3; sh makelibw3.sh'
    sys.stdout.write('executing "%s"\n' % strg)
    subprocess.call(strg,shell=True)
if not os.path.isfile('src/libwbacio_4.a'):
    strg = 'cd src/bacio; sh makebacio.sh'
    sys.stdout.write('executing "%s"\n' % strg)
    subprocess.call(strg,shell=True)

shtns_libdir = os.environ.get('SHTNS_LIBDIR')
fftw3_libdir = os.environ.get('FFTW3_LIBDIR')
if shtns_libdir is None or fftw3_libdir is None:
    raise ValueError('SHTNS_LIBDIR and FFTW3_LIBDIR env vars must be set')

srcs_read_sigma_spec =\
['src/_read_sigma_spec.pyf','src/kinds.f90','src/sigio_module.f90','src/shtns.f90','src/read_sigma_spec.f90']

libdirs = [shtns_libdir,fftw3_libdir]
ext_spec = Extension(name     = '_read_sigma_spec',
                     sources       = srcs_read_sigma_spec,
                     libraries     = ['shtns','fftw3'],
                     library_dirs  = libdirs)

srcs_read_sfc =\
['src/_read_sfc.pyf','src/kinds.f90','src/sfcio_module.f90','src/read_sfc.f90']

ext_sfc = Extension(name     = '_read_sfc',
                    sources       = srcs_read_sfc)

srcs_read_sigma_nemsio =\
['src/_read_sigma_nemsio.pyf','src/kinds.f90','src/nemsio_openclose.f90','src/nemsio_read.f90','src/nemsio_module.f90','src/read_sigma_nemsio.f90']

srcs_read_sfcflx_nemsio =\
['src/_read_sfcflx_nemsio.pyf','src/kinds.f90','src/nemsio_openclose.f90','src/nemsio_read.f90','src/nemsio_module.f90','src/read_sfcflx_nemsio.f90']

ext_sigma_nemsio = Extension(name     = '_read_sigma_nemsio',
                             sources       = srcs_read_sigma_nemsio,
                             libraries     = ['bacio_4','w3nco_d'],
                             library_dirs  = ['src'])

ext_sfcflx_nemsio = Extension(name     = '_read_sfcflx_nemsio',
                             sources       = srcs_read_sfcflx_nemsio,
                             libraries     = ['bacio_4','w3nco_d'],
                             library_dirs  = ['src'])

if __name__ == "__main__":
    setup(name = 'gfstonc',
          version           = "0.0.1",
          description       = "Modules and utilities for reading GFS output",
          author            = "Jeff Whitaker",
          author_email      = "jeffrey.s.whitaker@noaa.gov",
          url               = "http://github.com/jswhit/gfstonc",
          ext_modules       = [ext_spec,ext_sigma_nemsio,ext_sfcflx_nemsio,ext_sfc],
          packages          = ['ncepsigma','ncepsfc','ncepnemsio'],
          scripts           = ['utils/gfs_spectonc','utils/gfs_sfctonc','utils/gfs_nemsiotonc_3d','utils/gfs_nemsiotonc_2d'],
          )
