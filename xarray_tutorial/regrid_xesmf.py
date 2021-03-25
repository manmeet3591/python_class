ds_out = xr.Dataset({'lat': (['lat'], ds_anthrop.lat.values),
                     'lon': (['lon'], ds_anthrop.lon.values),
                    }
                   )
regridder = xe.Regridder(ds_imd_, ds_out, 'bilinear')
#regridder.clean_weight_file()
ds_imd_o = regridder(ds_imd_)
