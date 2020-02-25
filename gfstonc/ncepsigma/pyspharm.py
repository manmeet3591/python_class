import shtns
import numpy as np

class Spharmt(object):
    """
    wrapper class for commonly used spectral transform operations in
    atmospheric models.  Provides an interface to shtns compatible
    with pyspharm (pyspharm.googlecode.com).

    Jeffrey S. Whitaker <jeffrey.s.whitaker@noaa.gov>
    """
    def __init__(self,nlons,nlats,ntrunc,rsphere,gridtype='gaussian'):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes"""
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, \
                                shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
        if gridtype == 'gaussian':
            #self._shtns.set_grid(nlats,nlons,shtns.sht_gauss_fly|shtns.SHT_PHI_CONTIGUOUS,1.e-10)
            self._shtns.set_grid(nlats,nlons,shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,1.e-10)
        elif gridtype == 'regular':
            self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,1.e-10)
        self.lats = np.arcsin(self._shtns.cos_theta)
        self.lons = (2.*np.pi/nlons)*np.arange(nlons)
        self.nlons = nlons
        self.nlats = nlats
        self.ntrunc = ntrunc
        self.nlm = self._shtns.nlm
        self.degree = self._shtns.l
        self.lap = -self.degree*(self.degree+1.0).astype(np.complex)
        self.invlap = np.zeros(self.lap.shape, self.lap.dtype)
        self.invlap[1:] = 1./self.lap[1:]
        self.rsphere = rsphere
        self.lap = self.lap/rsphere**2
        self.invlap = self.invlap*rsphere**2
    def grdtospec(self,data):
        """compute spectral coefficients from gridded data"""
        data = np.ascontiguousarray(data, dtype=np.float)
        dataspec = np.empty(self.nlm, dtype=np.complex)
        self._shtns.spat_to_SH(data, dataspec)
        return dataspec
    def spectogrd(self,dataspec):
        """compute gridded data from spectral coefficients"""
        dataspec = np.ascontiguousarray(dataspec, dtype=np.complex)
        data = np.empty((self.nlats,self.nlons), dtype=np.float)
        self._shtns.SH_to_spat(dataspec, data)
        return data
    def getuv(self,vrtspec,divspec):
        """compute wind vector from spectral coeffs of vorticity and divergence"""
        vrtspec = np.ascontiguousarray(vrtspec, dtype=np.complex)
        divspec = np.ascontiguousarray(divspec, dtype=np.complex)
        u = np.empty((self.nlats,self.nlons), dtype=np.float)
        v = np.empty((self.nlats,self.nlons), dtype=np.float)
        self._shtns.SHsphtor_to_spat((self.invlap/self.rsphere)*vrtspec,\
               (self.invlap/self.rsphere)*divspec, u, v)
        return u,v
    def getvrtdivspec(self,u,v):
        """compute spectral coeffs of vorticity and divergence from wind vector"""
        u = np.ascontiguousarray(u, dtype=np.float)
        v = np.ascontiguousarray(v, dtype=np.float)
        vrtspec = np.empty(self.nlm, dtype=np.complex)
        divspec = np.empty(self.nlm, dtype=np.complex)
        self._shtns.spat_to_SHsphtor(u, v, vrtspec, divspec)
        return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec
    def getgrad(self,divspec):
        """compute gradient vector from spectral coeffs"""
        divspec = np.ascontiguousarray(divspec, dtype=np.complex)
        vrtspec = np.zeros(divspec.shape, dtype=np.complex)
        u = np.empty((self.nlats,self.nlons), dtype=np.float)
        v = np.empty((self.nlats,self.nlons), dtype=np.float)
        self._shtns.SHsphtor_to_spat(vrtspec,divspec)
        return u/rsphere,v/rsphere
