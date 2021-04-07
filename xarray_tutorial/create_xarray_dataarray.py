import numpy as np
import xarray as xr
import pandas as pd
times = pd.date_range("2000-01-01", periods=1)
data = np.random.rand(1, 4, 3)
lats = np.arange(0,4,1)
lons = np.arange(0,3,1)
foo = xr.DataArray(data, coords=[times, lats, lons], dims=["time", "lat", "lon"])
ds = xr.Dataset({
    'ndvi': xr.DataArray(
                data   = data,   # enter data here
                dims   = ['time'],
                coords = {'time': times},
                attrs  = {
                    '_FillValue': -999.9,
                    'units'     : 'W/m2'
                    }
                ),
            attrs = {'example_attr': 'this is a global attribute'}
    )
