import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import sys
from datetime import datetime
import xarray as xr
import pandas as pd
# Open the IMERG data file; it may be necessary to add a path if the file is not in the working directory.

file_ = sys.argv[1]
#f = h5py.File('3B-MO.MS.MRG.3IMERG.20140101-S000000-E235959.01.V06B.HDF5', 'r')
f = h5py.File(file_, 'r')

# View the available groups in the file and the variables in the 'Grid' group:

groups = [ x for x in f.keys() ]
print(groups)
gridMembers = [ x for x in f['Grid'] ]
print(gridMembers)

#print(file_[20:28],file_[31:35] )

precip = f['Grid/precipitationCal'][0][:][:]
precip = np.transpose(precip)
theLats = f['Grid/lat'][:]
theLons = f['Grid/lon'][:]
#theTime = f['Grid/time'][:]
#date_time_str = file_[26:30]+'-'+file_[30:32]+'-'+file_[32:34]+' '+file_[36:38]+'-'+file_[38:40]
#date_time_obj = datetime.strptime(date_time_str, '%Y-%m-%d %H-%M')
date_time_str = file_[26:30]+'-'+file_[30:32]+'-'+file_[32:34]+' '+file_[36:38]+':'+file_[38:40]
theTime =pd.Timestamp(date_time_str)
precip[precip==-9999.9]=np.nan
print(theTime)
ds = xr.Dataset({
    'pr': xr.DataArray(
        data   = precip[np.newaxis,:,:] ,   # enter data here
                dims   = ['time', 'lat', 'lon'],
                coords = {'time': np.atleast_1d(theTime), 'lat':theLats, 'lon':theLons},
                )
            },
    )
#print(ds)
ds.to_netcdf(file_[:-4]+'.nc')
#print(theLats, theLons, theTime, precip)
#print(theLats.shape, theLons.shape, theTime.shape, precip.shape)
#print(file_[26:34],file_[36:40] )
#print(date_time_obj)
