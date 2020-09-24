#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import stats
from datetime import datetime, timedelta
import cartopy.crs as ccrs
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
GeoAxes._pcolormesh_patched = Axes.pcolormesh
# fname = '/lus/dal/cccr_rnd/manmeet/iccp/data/INDIA.shp'
# adm1_shapes = list(shpreader.Reader(fname).geometries())
import xesmf as xe
import sys

data_dir_1 = '/home/cccr/msingh/dev_lab/ndvi3g_geo_v1_1981_0712.nc4'  # daily
data_dir_1 = sys.argv[1]  # daily
ds_gimms = xr.open_dataset(data_dir_1)

data_dir_2 = '/home/cccr/msingh/dev_lab/precip.1981.nc'
ds_cpc = xr.open_dataset(data_dir_2)


# In[2]:


# NDVI is a 15 day data and CPC is a monthly data
# Tasks to be performed: 
# 1. Regrid NDVI data to CPC
# 2. Negative to nan in NDVI, add time information in ndvi
# 3. CPC data to 15 day mean
# 4. Create event based time series in both the datasets
# 5. Use R code to do event coincidence analysis at 15 to 90 day coincidence period lag
# 6. Plot the Trigger coincidence rate and Precursor coincidence rate


# In[3]:


# But we are doing the analysis on monthly data following https://www.nature.com/articles/s41598-020-57910-1#MOESM1


# In[4]:


ndvi_ = ds_gimms.ndvi.values
ndvi_[ndvi_<0.0] = np.nan
ds_gimms['ndvi'] = (('time', 'lat', 'lon'),ndvi_)


# In[5]:


ds_out = xr.Dataset({'lat': (['lat'], ds_cpc.precip.lat.values),
                     'lon': (['lon'], ds_cpc.precip.lon.values),
                    }
                   )


# In[6]:


regridder = xe.Regridder(ds_gimms.ndvi, ds_out, 'bilinear')


# In[7]:


regridder.clean_weight_file()


# In[8]:


ds_ndvi_o = regridder(ds_gimms.ndvi)


# In[9]:


year = int(data_dir_1[-13:-9])


# In[25]:


import datetime
num=0
if data_dir_1[-9:-4]=='_0712':
    num=6
dt = [datetime.datetime(year, 1+num, 1, 0, 0),               datetime.datetime(year, 1+num, 15, 0, 0),               datetime.datetime(year, 2+num, 1, 0, 0),               datetime.datetime(year, 2+num, 15, 0, 0),               datetime.datetime(year, 3+num, 1, 0, 0),               datetime.datetime(year, 3+num, 15, 0, 0),               datetime.datetime(year, 4+num, 1, 0, 0),               datetime.datetime(year, 4+num, 15, 0, 0),               datetime.datetime(year, 5+num, 1, 0, 0),               datetime.datetime(year, 5+num, 15, 0, 0),               datetime.datetime(year, 6+num, 1, 0, 0),               datetime.datetime(year, 6+num, 15, 0, 0)]


# In[26]:


ds_gimms['time'] = ('time', dt)


# In[27]:


ds_gimms.to_netcdf(data_dir_1[-27:-4]+'_.nc', format="NETCDF3_CLASSIC")


# In[13]:


# save each file and then concatenate them, then use cdo monmax to use the 
# maximum value composite (MVC) technique as mentioned in the paper 
# https://www.nature.com/articles/s41598-020-57910-1#MOESM1 Chen et al 2020


# In[21]:


#data_dir_1[-27:-4]


# In[24]:


#data_dir_1[-9:-4]=='_0712'


# In[ ]:




