import os
import numpy as np
import xarray as xr

from datetime import timedelta, date



def daterange(start_date, end_date):
  for n in range(int((end_date - start_date).days)):
    yield start_date + timedelta(n)

start_date = date(2000, 2, 24)
#end_date = date(2020, 12, 31)
end_date = date(2000, 12, 31)
for single_date in daterange(start_date, end_date):
  print(single_date.strftime("%Y-%m-%d"), (single_date + timedelta(days=1)).strftime("%Y-%m-%d"))
  dates = single_date.strftime("%Y-%m-%d")
  datee = (single_date + timedelta(days=1)).strftime("%Y-%m-%d")
  os.system('python generate_modis_ndvi_gee.py '+dates+' '+datee+' &')
