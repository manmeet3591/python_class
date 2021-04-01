import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import stats
from datetime import datetime, timedelta
import cartopy.crs as ccrs
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.io.shapereader as shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
GeoAxes._pcolormesh_patched = Axes.pcolormesh
# fname = '/home/cccr/msingh/phd_obj2_data/ensemble_data/Admin2.shp'
# adm1_shapes = list(shapereader.Reader(fname).geometries())
# %matplotlib inline
import warnings
warnings.filterwarnings("ignore")
import xesmf as xe
from sklearn.metrics import mean_squared_error
import os
import matplotlib.pyplot as plt
import numpy.ma as ma
import seaborn as sns
import pandas as pd
import glob
import matplotlib.pyplot as plt
import skill_metrics as sm

data_dir = '/lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/CMIP6/CMIP6/'
ds_ipsl = xr.open_dataset(data_dir+'pr_day_IPSL-CM6A-LR_historical_r11i1p1f1_gr_18500101-20141231.nc').pr.resample(time="1MS").mean(dim="time")

ds_miroc6 = xr.open_mfdataset(data_dir+'pr_day_MIROC6_historical_r1i1p1f1_gn_*.nc', concat_dim='time')\
                        .pr.resample(time="1MS").mean(dim="time")

ds_noresm2 = xr.open_mfdataset(data_dir+'pr_day_NorESM2-LM_historical_r1i1p1f1_gn_*.nc', concat_dim='time')\
                        .pr.resample(time="1MS").mean(dim="time")

ds_ukesm1 = xr.open_mfdataset(data_dir+'pr_day_UKESM1-0-LL_historical_r1i1p1f2_gn_1*.nc', concat_dim='time')\
                        .pr.resample(time="1MS").mean(dim="time")

ds_mpiesm = xr.open_mfdataset(data_dir+'pr_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_*.nc', concat_dim='time')\
                        .pr.resample(time="1MS").mean(dim="time")

ds_cnrm = xr.open_mfdataset(data_dir+'pr_day_CNRM-ESM2-1_historical_r1i1p1f2_gr_1*.nc', concat_dim='time')\
                        .pr.resample(time="1MS").mean(dim="time")

data_dir = '/lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/IMD/imd-dly-p25-1901-2019.nc'
ds_imd = xr.open_dataset(data_dir).rf.resample(time="1MS").mean(dim="time")

lats, late, lons, lone = 14, 28, 74, 87
months = [6,7,8,9]

ds_ipsl_jjas_mt_ = ds_ipsl.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_ipsl.time.dt.year.isin(range(1901, 2015)))#.mean(dim='lat').mean(dim='lon')

ds_ipsl_jjas_mt = ds_ipsl_jjas_mt_.sel(time=ds_ipsl_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400

ds_miroc6_jjas_mt_ = ds_miroc6.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_miroc6.time.dt.year.isin(range(1901, 2015))).mean(dim='lat').mean(dim='lon')

ds_miroc6_jjas_mt = ds_miroc6_jjas_mt_.sel(time=ds_miroc6_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400

ds_noresm2_jjas_mt_ = ds_noresm2.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_noresm2.time.dt.year.isin(range(1901, 2015))).mean(dim='lat').mean(dim='lon')

ds_noresm2_jjas_mt = ds_noresm2_jjas_mt_.sel(time=ds_noresm2_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400

ds_ukesm1_jjas_mt_ = ds_ukesm1.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_ukesm1.time.dt.year.isin(range(1901, 2015))).mean(dim='lat').mean(dim='lon')

ds_ukesm1_jjas_mt = ds_ukesm1_jjas_mt_.sel(time=ds_ukesm1_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400

ds_mpiesm_jjas_mt_ = ds_mpiesm.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_mpiesm.time.dt.year.isin(range(1901, 2015))).mean(dim='lat').mean(dim='lon')

ds_mpiesm_jjas_mt = ds_mpiesm_jjas_mt_.sel(time=ds_mpiesm_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400

ds_cnrm_jjas_mt_ = ds_cnrm.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_cnrm.time.dt.year.isin(range(1901, 2015))).mean(dim='lat').mean(dim='lon')

ds_cnrm_jjas_mt = ds_cnrm_jjas_mt_.sel(time=ds_cnrm_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015))).mean(dim='lat').mean(dim='lon')

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time').values*86400


year = np.arange(1901,2015)

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015)))

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')


ds_ipsl_jjas_mt = ds_ipsl_jjas_mt_.sel(time=ds_ipsl_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_out = xr.Dataset({'lat': (['lat'], ds_ipsl_jjas_mt.lat.values),
                     'lon': (['lon'], ds_ipsl_jjas_mt.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_jjas_mt, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_jjas_mt)

print(np.corrcoef(ds_imd_o.mean(dim='year').values.flatten(), ds_ipsl_jjas_mt.mean(dim='year').values.flatten()),\
      np.std(ds_ipsl_jjas_mt.mean(dim='year').values.flatten())*86400,\
     centered_rms_dev(ds_ipsl_jjas_mt.mean(dim='year').values.flatten(),ds_imd_o.mean(dim='year').values.flatten()))

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015)))

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_miroc6_jjas_mt_ = ds_miroc6.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_miroc6.time.dt.year.isin(range(1901, 2015)))


ds_miroc6_jjas_mt = ds_miroc6_jjas_mt_.sel(time=ds_miroc6_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_out = xr.Dataset({'lat': (['lat'], ds_miroc6_jjas_mt.lat.values),
                     'lon': (['lon'], ds_miroc6_jjas_mt.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_jjas_mt, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_jjas_mt)

print(np.corrcoef(ds_imd_o.mean(dim='year').values.flatten(), ds_miroc6_jjas_mt.mean(dim='year').values.flatten()),\
      np.std(ds_miroc6_jjas_mt.mean(dim='year').values.flatten())*86400,\
     centered_rms_dev(ds_miroc6_jjas_mt.mean(dim='year').values.flatten(),ds_imd_o.mean(dim='year').values.flatten()))

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015)))

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_noresm2_jjas_mt_ = ds_noresm2.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_noresm2.time.dt.year.isin(range(1901, 2015)))


ds_noresm2_jjas_mt = ds_noresm2_jjas_mt_.sel(time=ds_noresm2_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_out = xr.Dataset({'lat': (['lat'], ds_noresm2_jjas_mt.lat.values),
                     'lon': (['lon'], ds_noresm2_jjas_mt.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_jjas_mt, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_jjas_mt)

print(np.corrcoef(ds_imd_o.mean(dim='year').values.flatten(), ds_noresm2_jjas_mt.mean(dim='year').values.flatten()),\
      np.std(ds_noresm2_jjas_mt.mean(dim='year').values.flatten())*86400,\
     centered_rms_dev(ds_noresm2_jjas_mt.mean(dim='year').values.flatten(),ds_imd_o.mean(dim='year').values.flatten()))

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015)))

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_ukesm1_jjas_mt_ = ds_ukesm1.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_ukesm1.time.dt.year.isin(range(1901, 2015)))


ds_ukesm1_jjas_mt = ds_ukesm1_jjas_mt_.sel(time=ds_ukesm1_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_out = xr.Dataset({'lat': (['lat'], ds_ukesm1_jjas_mt.lat.values),
                     'lon': (['lon'], ds_ukesm1_jjas_mt.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_jjas_mt, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_jjas_mt)

print(np.corrcoef(ds_imd_o.mean(dim='year').values.flatten(), ds_ukesm1_jjas_mt.mean(dim='year').values.flatten()),\
      np.std(ds_ukesm1_jjas_mt.mean(dim='year').values.flatten())*86400,\
     centered_rms_dev(ds_ukesm1_jjas_mt.mean(dim='year').values.flatten(),ds_imd_o.mean(dim='year').values.flatten()))

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015)))

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_mpiesm_jjas_mt_ = ds_mpiesm.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_mpiesm.time.dt.year.isin(range(1901, 2015)))


ds_mpiesm_jjas_mt = ds_mpiesm_jjas_mt_.sel(time=ds_mpiesm_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_out = xr.Dataset({'lat': (['lat'], ds_mpiesm_jjas_mt.lat.values),
                     'lon': (['lon'], ds_mpiesm_jjas_mt.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_jjas_mt, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_jjas_mt)

print(np.corrcoef(ds_imd_o.mean(dim='year').values.flatten(), ds_mpiesm_jjas_mt.mean(dim='year').values.flatten()),\
      np.std(ds_mpiesm_jjas_mt.mean(dim='year').values.flatten())*86400,\
     centered_rms_dev(ds_mpiesm_jjas_mt.mean(dim='year').values.flatten(),ds_imd_o.mean(dim='year').values.flatten()))

ds_imd_jjas_mt_ = ds_imd.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_imd.time.dt.year.isin(range(1901, 2015)))

ds_imd_jjas_mt = ds_imd_jjas_mt_.sel(time=ds_imd_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_cnrm_jjas_mt_ = ds_cnrm.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).\
sel(time=ds_cnrm.time.dt.year.isin(range(1901, 2015)))


ds_cnrm_jjas_mt = ds_cnrm_jjas_mt_.sel(time=ds_cnrm_jjas_mt_.time.dt.month.isin(months))\
.groupby('time.year').sum(dim='time')

ds_out = xr.Dataset({'lat': (['lat'], ds_cnrm_jjas_mt.lat.values),
                     'lon': (['lon'], ds_cnrm_jjas_mt.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_jjas_mt, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_jjas_mt)

print(np.corrcoef(ds_imd_o.mean(dim='year').values.flatten(), ds_cnrm_jjas_mt.mean(dim='year').values.flatten()),\
      np.std(ds_cnrm_jjas_mt.mean(dim='year').values.flatten())*86400,\
     centered_rms_dev(ds_cnrm_jjas_mt.mean(dim='year').values.flatten(),ds_imd_o.mean(dim='year').values.flatten()))

print(np.std(ds_imd_o.mean(dim='year').values.flatten()))

ccoef = np.array([1.0, 0.14, 0.29, -0.14, 0.35, -0.3, 0.43])
sdev = np.array([14.78, 9.60, 14.87, 9.87, 7.86, 13.43, 16.09])
crmsd = np.array([0.0, 12.62, 14.78,13.55,  14.99, 13.84, 14.78, ])

    # Specify labels for points in a cell array (M1 for model prediction 1,
    # etc.). Note that a label needs to be specified for the reference even
    # though it is not used.
label = ['Obs', 'IPSL', 'MIROC6', 'NORESM2', 'UKESM1', 'MPIESM', 'CNRM']

plt.figure(figsize=(12,8))
sm.taylor_diagram(sdev,crmsd,ccoef, markerLabel = label,
                      markerLabelColor = 'r', 
                      markerColor = 'r', markerLegend = 'on', 
                      tickRMS = range(0,20,4), 
                      colRMS = 'm', styleRMS = ':', widthRMS = 2.0, 
                      titleRMS = 'on', titleRMSDangle = 40.0, tickSTD = range(0,20,4),
                      axismax = 20.0, colSTD = 'b', styleSTD = '-.', 
                      widthSTD = 1.0, titleSTD = 'on', 
                      colCOR = 'k', styleCOR = '--', widthCOR = 1.0, 
                      titleCOR = 'on', markerSize = 10, alpha = 0.0)
