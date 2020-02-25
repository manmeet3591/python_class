import numpy as np
from .pyspharm import Spharmt
class ncepsigma(object):
    # read ncep 'sigma' file (spectral binary data)
    def __init__(self,filename):
        from _read_sigma_spec import read_specdata, read_header, read_griddata, get_vcoord
        nlons,nlats,nlevs,ntrunc,nvcoord,idate,fhour = read_header(filename)
        self.vcoord = get_vcoord(filename,nlevs,nvcoord).T
        self._read_specdata = read_specdata
        self._read_griddata = read_griddata
        self.nlons = nlons; self.nlats = nlats
        self.ntrunc = ntrunc; self.nlevs = nlevs
        self.idate = '%04i%02i%02i%02i' % (idate[3],idate[1],idate[2],idate[0])
        self.fhour = fhour
        self.filename = filename
        self.sp = Spharmt(nlons,nlats,ntrunc,6.3712e6,gridtype='gaussian')
        self._nf = np.sqrt(2.*np.pi)
        self.lats = (180./np.pi)*self.sp.lats
        self.lons = (360./nlons)*np.arange(nlons)
    def spectogrd(self,specdata):
        return self.sp.spectogrd(specdata)
    def getuv(self,vrtdata,divdata):
        return self.sp.getuv(vrtdata,divdata)
    def specdata(self):
        vrtspec, divspec,tempspec,zspec,lnpsspec,qspec,ozspec,cwmrspec =\
        self._read_specdata(self.filename,self.ntrunc,self.nlevs)
        nf = self._nf
        return nf*vrtspec.T,nf*divspec.T,nf*tempspec.T,\
               nf*zspec,nf*lnpsspec,nf*qspec.T,nf*ozspec.T,\
               nf*cwmrspec.T
    def griddata(self):
        ug,vg,tempg,zsg,psg,qg,ozg,cwmrg = self._read_griddata(self.filename,self.nlons,self.nlats,self.nlevs)
        return ug.T,vg.T,tempg.T,zsg.T,psg.T,qg.T,ozg.T,cwmrg.T
