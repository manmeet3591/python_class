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

lats, late, lons, lone = 14, 28, 74, 87
ds_anthrop_pr = []
ds_no_anthrop_pr = []
for i_y,year_ in enumerate(range(1982,2014)):
    print(year_)
####################################################################################
    data_dir = '/home/cccr/msingh/phd_data_cmip_account/anthrop/ATM/DAY/'+str(year_)
    filenames_a = glob.glob(data_dir+'/atm_*_day.nc')
    
    if len(filenames_a)==0:
            continue
    ens_list = []
    for file in filenames_a:
    
        #print(file)
        
        data_ = xr.open_dataset(file).pr.sel(lat=slice(late,lats))\
                        .sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon')
        datetimeindex = data_.indexes['time'].to_datetimeindex()
        data_['time'] = datetimeindex
        if data_.sel(time=slice(str(year_)+'-05-20', str(year_)+'-10-15')).time.shape[0]==148:
            ens_list.append(data_.sel(time=slice(str(year_)+'-05-20', str(year_)+'-10-15')))
    #print('test')
    if len(filenames_a)>0:
        ds_anthrop_pr_ = xr.concat(ens_list, dim='ens').mean(dim='ens')
        ds_anthrop_pr.append(ds_anthrop_pr_)
######################################################################################
    data_dir = '/home/cccr/msingh/phd_data_cmip_account/no_anthrop/ATM/DAY/'+str(year_)
    filenames_a = glob.glob(data_dir+'/atm_*_day.nc')
    if len(filenames_a)==0:
            continue
    ens_list = []
    #print(filenames_a)
    for file in filenames_a:
        #ens_list = []
        #print(file)
        data_ = xr.open_dataset(file).pr.sel(lat=slice(late,lats))\
                        .sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon')
        datetimeindex = data_.indexes['time'].to_datetimeindex()
        data_['time'] = datetimeindex
        if data_.sel(time=slice(str(year_)+'-05-20', str(year_)+'-10-15')).time.shape[0]==148:
            ens_list.append(data_.sel(time=slice(str(year_)+'-05-20', str(year_)+'-10-15')))
    if len(filenames_a)>0:
        ds_no_anthrop_pr_ = xr.concat(ens_list, dim='ens').mean(dim='ens')
        ds_no_anthrop_pr.append(ds_no_anthrop_pr_)
######################################################################################
ds_anthrop = xr.concat(ds_anthrop_pr, dim='time')
ds_no_anthrop = xr.concat(ds_no_anthrop_pr, dim='time')

data_dir = '/home/cccr/msingh/data/IMD/imd-dly-p25-1901-2019.nc'
ds_imd = xr.open_dataset(data_dir)
ds_imd['time'] = ds_imd.indexes['time'].normalize()
ds_anthrop['time'] = ds_anthrop.indexes['time'].normalize()
ds_no_anthrop['time'] = ds_no_anthrop.indexes['time'].normalize()
ds_imd_ = ds_imd.sel(time=ds_anthrop.time).sel(lat=slice(lats,late)).\
        sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon')
ds_anthrop_ = ds_anthrop*86400
ds_no_anthrop_ = ds_no_anthrop*86400

imd_active = np.zeros((2014-1982))
anthrop_active = np.zeros((2014-1982))
no_anthrop_active = np.zeros((2014-1982))

imd_break = np.zeros((2014-1982))
anthrop_break = np.zeros((2014-1982))
no_anthrop_break= np.zeros((2014-1982))
years = np.arange(1982,2014)
months=[6,7,8,9]
print(years.shape)

ds_imd_active= ds_imd_.rf.where(ds_imd_.rf>ds_imd_.rf.mean(dim='time')+ds_imd_.rf.std(dim='time'), drop=True)
ds_imd_active_ja = ds_imd_active.time.sel(time=ds_imd_active.time.dt.month.isin(months))
#print(ds_imd_active_ja.time)
cnt=0
for year in range(1982,2014):
    #print(ds_imd_active_ja.time.sel(time=ds_imd_active_ja.time.dt.year.isin([year])).values)
    #print(cnt)
    imd_active[cnt] = ds_imd_active_ja.time.sel(time=ds_imd_active_ja.time.dt.year.isin([year])).values.shape[0]
    cnt=cnt+1
#np.savetxt('imd_active.txt',ds_imd_active.time.sel(time=ds_imd_active.time.dt.year.isin([1982])).values )

ds_anthrop_active= ds_anthrop_.where(ds_anthrop_>ds_anthrop_.mean(dim='time')+ds_anthrop_.std(dim='time'), drop=True)

ds_anthrop_active_ja = ds_anthrop_active.time.sel(time=ds_anthrop_active.time.dt.month.isin(months))
cnt=0
for year in range(1982,2014):
    #print(ds_anthrop_active_ja.time.sel(time=ds_anthrop_active_ja.time.dt.year.isin([year])).values)
    anthrop_active[cnt] = ds_anthrop_active_ja.time.sel(time=ds_anthrop_active_ja.time.dt.year.isin([year])).values.shape[0]
    cnt=cnt+1
#ds_anthrop_.mean(dim='time')

ds_no_anthrop_active= ds_anthrop_.where(ds_no_anthrop_>ds_no_anthrop_.mean(dim='time')+ds_no_anthrop_.std(dim='time'), drop=True)

ds_no_anthrop_active_ja = ds_no_anthrop_active.time.sel(time=ds_no_anthrop_active.time.dt.month.isin(months))
cnt=0
for year in range(1982,2014):
    #print(ds_no_anthrop_active_ja.time.sel(time=ds_no_anthrop_active_ja.time.dt.year.isin([year])).values)
    no_anthrop_active[cnt] = ds_no_anthrop_active_ja.time.sel(time=ds_no_anthrop_active_ja.time.dt.year.isin([year])).values.shape[0]
    cnt=cnt+1
#ds_anthrop_.mean(dim='time')

ds_imd_break= ds_imd_.rf.where(ds_imd_.rf<ds_imd_.rf.mean(dim='time')-ds_imd_.rf.std(dim='time'), drop=True)
ds_imd_break_ja = ds_imd_break.time.sel(time=ds_imd_break.time.dt.month.isin(months))
#print(ds_imd_active_ja.time)
cnt=0
for year in range(1982,2014):
    #print(ds_imd_break_ja.time.sel(time=ds_imd_break_ja.time.dt.year.isin([year])).values)
    imd_break[cnt] = ds_imd_break_ja.time.sel(time=ds_imd_break_ja.time.dt.year.isin([year])).values.shape[0]
    cnt=cnt+1
#np.savetxt('imd_active.txt',ds_imd_active.time.sel(time=ds_imd_active.time.dt.year.isin([1982])).values )

ds_anthrop_break= ds_anthrop_.where(ds_anthrop_<ds_anthrop_.mean(dim='time')-ds_anthrop_.std(dim='time'), drop=True)

ds_anthrop_break_ja = ds_anthrop_break.time.sel(time=ds_anthrop_break.time.dt.month.isin(months))
cnt=0
for year in range(1982,2014):
    #print(ds_anthrop_break_ja.time.sel(time=ds_anthrop_break_ja.time.dt.year.isin([year])).values)
    anthrop_break[cnt] = ds_anthrop_break_ja.time.sel(time=ds_anthrop_break_ja.time.dt.year.isin([year])).values.shape[0]
    cnt=cnt+1
#ds_anthrop_.mean(dim='time')

ds_no_anthrop_break= ds_no_anthrop_.where(ds_no_anthrop_<ds_no_anthrop_.mean(dim='time')-ds_no_anthrop_.std(dim='time'), drop=True)

ds_no_anthrop_break_ja = ds_no_anthrop_break.time.sel(time=ds_no_anthrop_break.time.dt.month.isin(months))
cnt=0
for year in range(1982,2014):
    #print(ds_no_anthrop_break_ja.time.sel(time=ds_no_anthrop_break_ja.time.dt.year.isin([year])).values)
    no_anthrop_break[cnt] = ds_no_anthrop_break_ja.time.sel(time=ds_no_anthrop_break_ja.time.dt.year.isin([year])).values.shape[0]
    cnt=cnt+1
#ds_no_anthrop_.mean(dim='time')

data = {'years':years, 'data': imd_active}
df = pd.DataFrame(data=data)
df['kind'] = 'IMD'

data = {'years':years, 'data': anthrop_active}
df1 = pd.DataFrame(data=data)
df1['kind'] = 'Anthrop'

data = {'years':years, 'data': no_anthrop_active}
df2 = pd.DataFrame(data=data)
df2['kind'] = 'No-anthrop'

df_active = df.append(df1).append(df2)

fig = plt.figure(figsize=(20,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
g = sns.barplot(x="years", y="data", hue="kind", data=df_active, ax=ax)
g.legend_.set_title(None)
plt.grid()
plt.ylabel('Number of active days in JJAS', fontsize=18)
plt.xlabel('Years', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14, rotation=45)
plt.legend(fontsize=18)

data = {'years':years, 'data': imd_break}
df = pd.DataFrame(data=data)
df['kind'] = 'IMD'

data = {'years':years, 'data': anthrop_break}
df1 = pd.DataFrame(data=data)
df1['kind'] = 'Anthrop'

data = {'years':years, 'data': no_anthrop_break}
df2 = pd.DataFrame(data=data)
df2['kind'] = 'No-anthrop'

df_break = df.append(df1).append(df2)

fig = plt.figure(figsize=(20,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
g = sns.barplot(x="years", y="data", hue="kind", data=df_break, ax=ax)
g.legend_.set_title(None)
plt.grid()
plt.ylabel('Number of break days in JJAS', fontsize=18)
plt.xlabel('Years', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14, rotation=45)
plt.legend(fontsize=18)
