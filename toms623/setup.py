from numpy.distutils.core  import setup, Extension

# interface for Renka's algorithm 623 fortran code
ext623 = Extension(name  = '_toms623',
                sources       = ['_toms623.pyf','_toms623.f90'])
ext772 = Extension(name  = '_stripack',
                sources       = ['_stripack.pyf','_stripack.f'])


if __name__ == "__main__":
    setup(name = 'toms623',
          version           = "0.9",
          description       = "Python interface to TOMS 623 fortran code",
          author            = "Jeff Whitaker",
          author_email      = "jeffrey.s.whitaker@noaa.gov",
          ext_modules       = [ext623,ext772],
          packages          = ['toms623'],
          )
