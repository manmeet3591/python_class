import _toms623, _stripack
import numpy as np

class trintrp(object):
    def __init__(self, lons, lats, reorder = None):
        """
given mesh points (lons,lats in radians) define triangulation.
n is size of input mesh (length of 1-d arrays lons and lats).
triangulation can sometimes be sped up be re-ordering mesh coordinates
using reorder kwarg (options are 'lon','lat', and None (default)).
reorder='lat' causes lon,lat pairs to be reordered by increasing
latitude before triangulation is computed.

Algorithm:
 R. J. Renka, "ALGORITHM 623:  Interpolation on the Surface of a
 Sphere", ACM Trans. Math. Software, Vol. 10, No. 4, December 1984,
 pp. 437-439."""
        if len(lons.shape) != 1 or len(lats.shape) != 1:
            raise ValueError('lons and lats must be 1d')
        lats = lats.astype(np.float64,copy=False)
        lons = lons.astype(np.float64,copy=False)
        if (np.abs(lons)).max() > 2.*np.pi:
            msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
            raise ValueError(msg)
        if (np.abs(lats)).max() > 0.5*np.pi:
            msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
            raise ValueError(msg)
        npts = len(lons)
        if len(lats) != npts:
            raise ValueError('lons and lats must have same length')
        # compute cartesian coords on unit sphere.
        x,y,z = _toms623.trans(lats,lons,npts)
        # reorder the coordinates by latitude.
        if reorder == 'lat':
            dummy = lats.copy()
        elif reorder == 'lon':
            dummy = lons.copy()
        elif reorder is None:
            dummy = None
        else:
            raise ValueError('illegal value for reorder kwarg')
        if dummy is not None:
            self.ind = _toms623.reordr(4,dummy,x,y,z,npts)
        else:
            self.ind = np.arange(1,npts+1,1)
        self.ind = self.ind - 1 # change to zero based indexing
        # create the triangulation.
        iadj,iend,ierr = _toms623.trmesh(x,y,z,npts)
        if ierr != 0:
            raise ValueError('ierr = %s in trmesh' % ierr)
        self.lons = lons; self.lats = lats; self.npts = npts
        self.x = x; self.y = y; self.z = z
        self.iadj = iadj; self.iend = iend
    def interp(self,olons,olats,data,order=1):
        """
given a triangulation, perform interpolation on
oints olons,olats (in radians), return result in data.
olons, olats can be 1d or 2d (output data array has same shape as olats,lons).
order of interpolation specified by 'order' kwarg, can be 0 (nearest neighbor),
1 (linear), or 3 (hermite cubic). Default is linear."""
        shapeout = olons.shape
        if len(shapeout) not in [1,2]:
            raise ValueError('olons,olats must be 1d or 2d')
        olons1 = (olons.astype(np.float64,copy=False)).ravel()
        olats1 = (olats.astype(np.float64,copy=False)).ravel()
        nptso = len(olons1)
        if (np.abs(olons1)).max() > 2.*np.pi:
            msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
            raise ValueError(msg)
        if (np.abs(olats1)).max() > 0.5*np.pi:
            msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
            raise ValueError(msg)
        if len(olats1) != nptso:
            raise ValueError('lons and lats must have same length')
        if len(data) != self.npts:
            raise ValueError('input data wrong size')
        # reorder input data based on sorting of nodes.
        data_reordered = data[self.ind].astype(np.float64,copy=False)
        if order == 0:
            odata,ierr = \
            _toms623.intrpnn_n(olats1, olons1,\
                         self.x, self.y, self.z, data_reordered,\
                         self.iadj,self.iend,self.npts,nptso)
        elif order == 1:
            odata,ierr = \
            _toms623.intrpc0_n(olats1, olons1,\
                         self.x, self.y, self.z, data_reordered,\
                         self.iadj,self.iend,self.npts,nptso)
        elif order == 3:
            odata,ierr = \
            _toms623.intrpc1_n(olats1, olons1,\
                         self.x, self.y, self.z, data_reordered,\
                         self.iadj,self.iend,self.npts,nptso)
        else:
            raise ValueError('order must be 0,1 or 3')
        if ierr != 0:
            raise ValueError('ierr = %s in intrpc0_n' % ierr)
        return odata.reshape(shapeout)
    def interp_nn(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=0)"""
        return self.interp(olons,olats,data,order=0)
    def interp_linear(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=1)"""
        return self.interp(olons,olats,data,order=1)
    def interp_cubic(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=3)"""
        return self.interp(olons,olats,data,order=3)

class stripack(object):
    def __init__(self, lons, lats):
        """
given mesh points (lons,lats in radians) define triangulation.
n is size of input mesh (length of 1-d arrays lons and lats).

Same as trintrp, but uses improved triangulation routines
from TOMS 772, and has no 'reorder' kwarg.

Algorithm:
 R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay triangulation
 and Voronoi diagram on the surface of a sphere"
 ACM Trans. Math. Software, Volume 23 Issue 3, Sept. 1997
 pp 416-434."""
        if len(lons.shape) != 1 or len(lats.shape) != 1:
            raise ValueError('lons and lats must be 1d')
        lats = lats.astype(np.float64,copy=False)
        lons = lons.astype(np.float64,copy=False)
        if (np.abs(lons)).max() > 2.*np.pi:
            msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
            raise ValueError(msg)
        if (np.abs(lats)).max() > 0.5*np.pi:
            msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
            raise ValueError(msg)
        npts = len(lons)
        if len(lats) != npts:
            raise ValueError('lons and lats must have same length')
        # compute cartesian coords on unit sphere.
        x,y,z = _stripack.trans(lats,lons,npts)
        lst,lptr,lend,ierr = _stripack.trmesh(x,y,z,npts)
        if ierr != 0:
            raise ValueError('ierr = %s in trmesh' % ierr)
        self.lons = lons; self.lats = lats; self.npts = npts
        self.x = x; self.y = y; self.z = z
        self.lptr = lptr; self.lst = lst; self.lend = lend
    def interp(self,olons,olats,data,order=1):
        """
given a triangulation, perform interpolation on
oints olons,olats (in radians), return result in data.
olons, olats can be 1d or 2d (output data array has same shape as olats,lons).
order of interpolation specified by 'order' kwarg, can be 0 (nearest neighbor),
1 (linear), or 3 (hermite cubic, no tension). Default is linear."""
        shapeout = olons.shape
        if len(shapeout) not in [1,2]:
            raise ValueError('olons,olats must be 1d or 2d')
        olons1 = (olons.astype(np.float64,copy=False)).ravel()
        olats1 = (olats.astype(np.float64,copy=False)).ravel()
        nptso = len(olons1)
        if (np.abs(olons1)).max() > 2.*np.pi:
            msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
            raise ValueError(msg)
        if (np.abs(olats1)).max() > 0.5*np.pi:
            msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
            raise ValueError(msg)
        if len(olats1) != nptso:
            raise ValueError('lons and lats must have same length')
        if len(data) != self.npts:
            raise ValueError('input data wrong size')
        if order == 0:
            odata,ierr = \
            _stripack.intrpnn_n(olats1, olons1,\
                         self.x, self.y, self.z, data.astype(np.float64),\
                         self.lst,self.lptr,self.lend,self.npts,nptso)
        elif order == 1:
            odata,ierr = \
            _stripack.intrpc0_n(olats1, olons1,\
                         self.x, self.y, self.z, data.astype(np.float64),\
                         self.lst,self.lptr,self.lend,self.npts,nptso)
        elif order == 3:
            odata,ierr = \
            _stripack.intrpc1_n(olats1, olons1,\
                         self.x, self.y, self.z, data.astype(np.float64),\
                         self.lst,self.lptr,self.lend,self.npts,nptso)
        else:
            raise ValueError('order must be 0,1 or 3')
        if ierr != 0:
            raise ValueError('ierr = %s in intrpc0_n' % ierr)
        return odata.reshape(shapeout)
    def interp_nn(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=0)"""
        return self.interp(olons,olats,data,order=0)
    def interp_linear(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=1)"""
        return self.interp(olons,olats,data,order=1)
    def interp_cubic(self,olons,olats,data):
        """
same as interp(olons,olats,data,order=3)"""
        return self.interp(olons,olats,data,order=3)
