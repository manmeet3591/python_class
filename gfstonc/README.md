# gfstonc
python utilities and modules to read NCEP gfs spectral, surface and nemsio files and convert to netcdf.

[numpy](http://numpy.org), [netcdf4-python](https://github.com/Unidata/netcdf4-python) and
[gfortran](https://gcc.gnu.org/wiki/GFortran) are required, 
plus the [shtns](https://bitbucket.org/nschaeff/shtns) spherical harmonic transform
library (including the python interface).

* set the env vars `SHTNS_LIBDIR` and `FFTW3_LIBDIR` to point to where the `shtns` and `fftw3` 
libraries are installed.  [fftw3](http://www.fftw.org) is required by [shtns](https://bitbucket.org/nschaeff/shtns).
* `python setup.py build`
   - setup.py will try to build `src/libw3nco_d.a`  and `src/libbacio_4.a` if they do not
already exist. 
* `python setup.py install` (or `python setup.py install --user` to install in your 
home directory).

*Will not work on Windows!*

`utils/gfs_spectonc` converts gfs binary spectral files to netcdf (assumes the binary spectral data is big endian).

`utils/gfs_sfctonc` converts gfs binary surface files to netcdf (assumes the binary surface data is big endian).

`utils/gfs_nemsiotonc_3d` converts gfs 3d nemsio files to netcdf.

`utils/gfs_nemsiotonc_2d` converts gfs 2d nemsio files to netcdf.

Watch where setup.py installs the utility scripts, you will need to add that location to your `$PATH`.
