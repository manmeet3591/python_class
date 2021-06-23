import h5py
import numpy as np
import datetime as dt
import xarray as xr
import sys

input_file = sys.argv[1]
f = h5py.File(input_file, 'r')

yy = np.max(f['Grid']['GridTime']['Year'][:])
mm = np.max(f['Grid']['GridTime']['Month'][:])
dd = np.max(f['Grid']['GridTime']['DayOfMonth'][:])

time_ = dt.datetime(yy, mm, dd)
lat = np.linspace(-67+0.25/2,67-0.25/2,536)
lon = np.linspace(-180+0.25/2,180-0.25/2,1440)
lev = np.arange(0.125,0.125+0.25*80,0.25)

lh = f['Grid']['latentHeating'][:]
spr = f['Grid']['surfacePrecipRate'][:]
stratfrac = f['Grid']['stratiformFraction'][:]
lh[lh==-9999.9] = np.nan
spr[spr==-9999.9] = np.nan
stratfrac[stratfrac==-9999.9] = np.nan

stratfrac = stratfrac.T
spr = spr.T
lh = lh.swapaxes(1,2)
lh = lh[np.newaxis,...]

stratfrac = stratfrac[np.newaxis,...]
spr = spr[np.newaxis,...]

df = xr.Dataset(
        data_vars={'latentHeating':    (('time', 'level', 'latitude', 'longitude' ), lh),
                   'surfacePrecipRate':    (('time', 'latitude', 'longitude' ), spr),
                   'stratiformFraction':    (('time', 'latitude', 'longitude'), stratfrac)},
        coords={ 'time':np.atleast_1d((time_-dt.datetime(1,1,1)).days+2),
                'level': lev,
                'longitude': lon,
                'latitude': lat,
            })
df.latentHeating.attrs['long_name'] = 'Latent Heating'
df.surfacePrecipRate.attrs['long_name'] = 'Surface Precipitation Rate'
df.stratiformFraction.attrs['long_name'] = 'Stratiform Fraction'

df.time.attrs['units'] = 'days since 01-01-01 00:00:00'
df.time.attrs['standard_name'] = 'time'
df.time.attrs['long_name'] = 'Year AD'
df.time.attrs['calendar'] = 'standard'

df.longitude.attrs['standard_name'] = 'longitude'
df.longitude.attrs['long_name'] = 'longitude'
df.longitude.attrs['units'] = 'degrees_east'
df.longitude.attrs['axis'] = 'X'
df.latitude.attrs['standard_name'] = 'latitude'
df.latitude.attrs['long_name'] = 'latitude'
df.latitude.attrs['units'] = 'degrees_north'
df.latitude.attrs['axis'] = 'Y'

df.level.attrs['standard_name'] = 'air_pressure'
df.level.attrs['long_name'] = 'pressure_level'
df.level.attrs['units'] = 'km'
df.level.attrs['axis'] = 'Z'
df.level.attrs['positive'] = 'up'

df.to_netcdf(input_file[:-4]+'nc')
