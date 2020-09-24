import numpy as np
class ncepnemsio_3d(object):
    # read 3d ncep gfs 'nemsio' file
    def __init__(self,filename):
        from _read_sigma_nemsio import read_nemsio_header, read_nemsio_griddata, read_nemsio_coords
        nlons,nlats,nlevs,idate,nfhour = read_nemsio_header(filename)
        vcoord,lats,lons = read_nemsio_coords(filename,nlons,nlats,nlevs)
        self.vcoord = vcoord[:,:2,0].T
        self._read_griddata = read_nemsio_griddata
        self.nlons = nlons; self.nlats = nlats
        self.nlevs = nlevs
        self.idate = '%04i%02i%02i%02i' % (idate[0],idate[1],idate[2],idate[3])
        self.fhour = nfhour
        self.filename = filename
        self.lats = lats
        self.lons = lons
    def griddata(self):
        ug,vg,tempg,zsg,psg,qg,ozg,cwmrg,dpresg,presg = self._read_griddata(self.filename,self.nlons,self.nlats,self.nlevs)
        return ug.T,vg.T,tempg.T,zsg.T,psg.T,qg.T,ozg.T,cwmrg.T,dpresg.T,presg.T
class ncepnemsio_2d(object):
    # read ncep gfs sfc or flx 'nemsio' file with only 2d fields.
    def __init__(self,filename):
        from _read_sfcflx_nemsio import read_nemsio_header,\
                                        read_nemsio_varnames,\
                                        read_nemsio_2dgriddata, read_nemsio_latlons
        nlons,nlats,nrecs,idate,nfhour = read_nemsio_header(filename)
        lats,lons = read_nemsio_latlons(filename,nlons,nlats)
        irecnames, ireclevtypes, ireclevs = read_nemsio_varnames(filename,nrecs)
        self._read_griddata = read_nemsio_2dgriddata
        self.nlons = nlons; self.nlats = nlats
        self.nrecs = nrecs
        self.idate = '%04i%02i%02i%02i' % (idate[0],idate[1],idate[2],idate[3])
        self.fhour = nfhour
        self.filename = filename
        self.lats = lats
        self.lons = lons
        self.irecnames = irecnames
        self.ireclevtypes = ireclevtypes
        self.reclevs = ireclevs
        recnames = []; reclevtypes = []
        for n in range(nrecs):
            s = irecnames[n].tostring()
            s = s.encode('ascii').replace('\x00','').strip()
            recnames.append(s)
            s = ireclevtypes[n].tostring()
            s = s.encode('ascii').replace('\x00','').strip()
            reclevtypes.append(s)
        self.recnames = recnames
        self.reclevtypes = reclevtypes
    def griddata(self):
        grids = self._read_griddata(self.filename,self.nlons,self.nlats,self.irecnames,self.ireclevtypes, self.reclevs)
        return grids.T
