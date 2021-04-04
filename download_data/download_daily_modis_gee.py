import ee
import numpy as np
ee.Initialize()
import  xarray as xr
#data_dir = '/home/cccr/msingh/data/merra_aod_global_1980_2020/MERRA2_100.tavg1_2d_aer_Nx.1981_monmean.nc4'
import sys
dates = sys.argv[1]# '2020-03-23'
datee = sys.argv[2]# '2020-04-24'
#ds_merra2 = xr.open_dataset(data_dir)
collection_evi = ee.ImageCollection('MODIS/MYD09GA_006_EVI').select('EVI').filterDate(dates,datee ).mean()
collection_ndsi = ee.ImageCollection('MODIS/MYD09GA_006_NDSI').select('NDSI').filterDate(dates,datee ).mean()
collection_ndwi = ee.ImageCollection('MODIS/MYD09GA_006_NDWI').select('NDWI').filterDate(dates,datee ).mean()
collection_ndvi = ee.ImageCollection('MODIS/MYD09GA_006_NDVI').select('NDVI').filterDate(dates,datee ).mean()
#p = ee.Geometry.Point(32.3, 40.3)
#lats = ds_merra2.lat.values
#lons = ds_merra2.lon.values
lats = np.arange(14,28,0.25)
lons = np.arange(74,87,0.25)
#print(lats.shape, lons.shape, ds_merra2.TOTEXTTAU.values.shape)
evi = np.zeros((lats.shape[0], lons.shape[0]))
ndsi = np.zeros((lats.shape[0], lons.shape[0]))
ndwi = np.zeros((lats.shape[0], lons.shape[0]))
ndvi = np.zeros((lats.shape[0], lons.shape[0]))

#data = collection_evi.reduceRegion(ee.Reducer.first(),p,50000).get("EVI")# 0.5 degree = 50km =50000
#dataN = ee.Number(data)
#print(type(dataN.getInfo()))
#print(dataN.getInfo())
for i_lat in range(lats.shape[0]):
  for j_lon in range(lons.shape[0]):
    p = ee.Geometry.Point(lons[j_lon], lats[i_lat]) #p = ee.Geometry.Point([lon, lat])
    data = collection_evi.reduceRegion(ee.Reducer.mean(),p,25000).get("EVI")# 0.5 degree = 50km =50000
    evi[i_lat, j_lon] = ee.Number(data).getInfo()
    data = collection_ndsi.reduceRegion(ee.Reducer.mean(),p,25000).get("NDSI")# 0.5 degree = 50km =50000
    ndsi[i_lat, j_lon] = ee.Number(data).getInfo()
    data = collection_ndwi.reduceRegion(ee.Reducer.mean(),p,25000).get("NDWI")# 0.5 degree = 50km =50000
    ndwi[i_lat, j_lon] = ee.Number(data).getInfo()
    data = collection_ndvi.reduceRegion(ee.Reducer.mean(),p,25000).get("NDVI")# 0.5 degree = 50km =50000
    ndvi[i_lat, j_lon] = ee.Number(data).getInfo()
    print(lats[i_lat], lons[j_lon], ndvi[i_lat, j_lon])
