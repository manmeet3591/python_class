lats_ = np.unique(np.asarray(lat_))
lons_ = np.unique(np.asarray(lon_))
dates_ = np.unique(np.asarray(date_))
data_np = np.empty((dates_.shape[0], lats_.shape[0], lons_.shape[0]))
ds_jf = xr.DataArray(data=data_np, dims=['time', 'lat', 'lon'], coords={'time':dates_, 'lat':lats_, 'lon':lons_})
ds_jm = xr.DataArray(data=data_np, dims=['time', 'lat', 'lon'], coords={'time':dates_, 'lat':lats_, 'lon':lons_})
ds_sam = xr.DataArray(data=data_np, dims=['time', 'lat', 'lon'], coords={'time':dates_, 'lat':lats_, 'lon':lons_})
ds_af = xr.DataArray(data=data_np, dims=['time', 'lat', 'lon'], coords={'time':dates_, 'lat':lats_, 'lon':lons_})
ds_am = xr.DataArray(data=data_np, dims=['time', 'lat', 'lon'], coords={'time':dates_, 'lat':lats_, 'lon':lons_})
for i_ in range(df_cclm_mpiesm.shape[0]):
    ds_jf.loc[dict(lat=lat_[i_], lon=lon_[i_], time=date_[i_])] = jf[i_]
    ds_jm.loc[dict(lat=lat_[i_], lon=lon_[i_], time=date_[i_])] = jm[i_]
    ds_sam.loc[dict(lat=lat_[i_], lon=lon_[i_], time=date_[i_])] = sam[i_]
    ds_af.loc[dict(lat=lat_[i_], lon=lon_[i_], time=date_[i_])] = af[i_]
    ds_am.loc[dict(lat=lat_[i_], lon=lon_[i_], time=date_[i_])] = am[i_]
