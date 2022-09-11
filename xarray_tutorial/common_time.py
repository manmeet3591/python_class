for arr in [gfs_f24_e,ds_tp_era5, ds_dl, ds_lr]:
    arr['time'] = arr.indexes['time'].normalize()

ds_tp_era5_e = ds_tp_era5.sel(time=ds_tp_era5.time.isin(gfs_f24_e.time.values))
ds_dl_e = ds_dl.sel(time=ds_dl.time.isin(gfs_f24_e.time.values))
ds_lr_e = ds_lr.sel(time=ds_lr.time.isin(gfs_f24_e.time.values))
