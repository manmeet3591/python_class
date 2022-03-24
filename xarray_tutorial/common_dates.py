ds_obs = ds_observed[i_var]

file_ = base_dir+'historical_v1/'+var+'_historical_'+df_corr_['model'].iloc[cnt]+'_ens_'+str(df_corr_['ensemble'].iloc[cnt])+'.nc'
print(file_)
ds_historical = xr.open_dataset(file_).sel(lat=location_[0], lon= 360-location_[1], method='nearest')
datetimeindex = ds_historical.indexes['time'].to_datetimeindex()
ds_historical['time'] = datetimeindex

common_dates = np.intersect1d(np.asarray(ds_historical.time.values.astype('datetime64[D]').tolist()), \
              np.asarray(ds_obs.time.values.astype('datetime64[D]').tolist()))
ds_historical['time'] = ('time', ds_historical.time.values.astype('datetime64[D]'))
ds_historical_ = ds_historical.sel(time=common_dates)
ds_obs_ = ds_obs.sel(time=common_dates)

file_ = base_dir+'ssp126_v1/'+var+'_ssp126_'+df_corr_['model'].iloc[cnt]+'_ens_'+str(df_corr_['ensemble'].iloc[cnt])+'.nc'
print(file_)
ds_ssp126 = xr.open_dataset(file_).sel(lat=location_[0], lon= 360-location_[1], method='nearest')

file_ = base_dir+'ssp585_v1/'+var+'_ssp585_'+df_corr_['model'].iloc[cnt]+'_ens_'+str(df_corr_['ensemble'].iloc[cnt])+'.nc'
print(file_)
ds_ssp585 = xr.open_dataset(file_).sel(lat=location_[0], lon= 360-location_[1], method='nearest')
