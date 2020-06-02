import numpy as np
from ncepsigma import Spharmt
class ncepsfc(object):
    # read ncep 'sfc' file (fortran gridded binary data)
    def __init__(self,filename):
        from _read_sfc import read_griddata, read_header
        nlons,nlats,lsoil,idate,fhour = read_header(filename)
        self._read_griddata = read_griddata
        self.nlons = nlons; self.nlats = nlats
        self.lsoil = lsoil
        self.idate = '%04i%02i%02i%02i' % (idate[3],idate[1],idate[2],idate[0])
        self.fhour = fhour
        self.filename = filename
        sp = Spharmt(nlons,nlats,nlats/2,6.3712e6,gridtype='gaussian')
        self.lats = (180./np.pi)*sp.lats
        self.lons = (360./nlons)*np.arange(nlons)
    def griddata(self):
        grids2d,grids2d_desc,grids2d_name,grids3d,grids3d_desc,grids3d_name = self._read_griddata(self.filename,self.nlons,self.nlats,self.lsoil)
        grds2d_desc = []
        for n in range(grids2d_desc.shape[0]):
            s = grids2d_desc[n].tostring()
            s = s.encode('ascii').replace('\x00','').strip()
            grds2d_desc.append(s)
        grds2d_name = []
        for n in range(grids2d_name.shape[0]):
            s = grids2d_name[n].tostring()
            s = s.encode('ascii').replace('\x00','').strip()
            grds2d_name.append(s)
        grds3d_desc = []
        for n in range(grids3d_desc.shape[0]):
            s = grids3d_desc[n].tostring()
            s = s.encode('ascii').replace('\x00','').strip()
            grds3d_desc.append(s)
        grds3d_name = []
        for n in range(grids3d_name.shape[0]):
            s = grids3d_name[n].tostring()
            s = s.encode('ascii').replace('\x00','').strip()
            grds3d_name.append(s)
        return grids2d.T,grds2d_desc,grds2d_name,grids3d.T,grds3d_desc,grds3d_name
