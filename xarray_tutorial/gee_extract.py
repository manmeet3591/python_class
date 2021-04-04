import ee
import numpy as np
ee.Initialize()
import  xarray as xr
import sys
#data_dir = '/home/cccr/msingh/data/merra_aod_global_1980_2020/MERRA2_100.tavg1_2d_aer_Nx.1981_monmean.nc4'
import sys
#dates = sys.argv[1]# '2020-03-23'
#datee = sys.argv[2]# '2020-04-24'
#ds_merra2 = xr.open_dataset(data_dir)
from datetime import timedelta, date
import pandas as pd


def daterange(start_date, end_date):
  for n in range(int((end_date - start_date).days)):
    yield start_date + timedelta(n)

start_date = date(2003, 1, 1)
end_date = date(2021, 1, 1)
#end_date = date(2000, 12, 31)
#end_date = date(2021, 1, 1)


times = pd.date_range('2003-01-01', '2021-01-01')
print(times, times.shape)
#sys.exit()
evi  = np.zeros((times.shape[0]))
ndsi = np.zeros((times.shape[0]))
ndwi = np.zeros((times.shape[0]))
ndvi = np.zeros((times.shape[0]))
monsoon_region = ee.Geometry.Polygon([

  [[74, 14], [87, 14],[87, 28], [74, 28], [74, 14]]

  ]);


cnt = 0
for single_date in daterange(start_date, end_date):
  print(single_date.strftime("%Y-%m-%d"), (single_date + timedelta(days=1)).strftime("%Y-%m-%d"))
  dates = single_date.strftime("%Y-%m-%d")
  datee = (single_date + timedelta(days=1)).strftime("%Y-%m-%d")


  collection_evi = ee.ImageCollection('MODIS/MYD09GA_006_EVI').select('EVI').filterDate(dates,datee ).mean()
  collection_ndsi = ee.ImageCollection('MODIS/MYD09GA_006_NDSI').select('NDSI').filterDate(dates,datee ).mean()
  collection_ndwi = ee.ImageCollection('MODIS/MYD09GA_006_NDWI').select('NDWI').filterDate(dates,datee ).mean()
  collection_ndvi = ee.ImageCollection('MODIS/MYD09GA_006_NDVI').select('NDVI').filterDate(dates,datee ).mean()

  data = collection_evi.reduceRegion(ee.Reducer.mean(), monsoon_region, scale=500).get('EVI')
  evi[cnt] = ee.Number(data).getInfo()

  data = collection_ndsi.reduceRegion(ee.Reducer.mean(), monsoon_region, scale=500).get('NDSI')
  ndsi[cnt] = ee.Number(data).getInfo()

  data = collection_ndwi.reduceRegion(ee.Reducer.mean(), monsoon_region, scale=500).get('NDWI')
  ndwi[cnt] = ee.Number(data).getInfo()

  data = collection_ndvi.reduceRegion(ee.Reducer.mean(), monsoon_region, scale=500).get('NDVI')
  ndvi[cnt] = ee.Number(data).getInfo()

  cnt = cnt+1
np.savetxt('evi.out', evi, delimiter=',')
np.savetxt('ndsi.out', ndsi, delimiter=',')
np.savetxt('ndwi.out', ndwi, delimiter=',')
np.savetxt('ndvi.out', ndvi, delimiter=',')
ds = xr.DataArray(evi, coords=[times], dims=["time"])
ds.to_netcdf('evi_.nc')
ds = xr.DataArray(ndsi, coords=[times], dims=["time"])
ds.to_netcdf('ndsi_.nc')
ds = xr.DataArray(ndwi, coords=[time], dims=["time"])
ds.to_netcdf('ndwi_.nc')
ds = xr.DataArray(ndvi, coords=[times], dims=["time"])
ds.to_netcdf('ndvi_.nc')
