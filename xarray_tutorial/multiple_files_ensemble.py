data_dir = '/home/cccr/msingh/phd_data_cmip_account/anthrop/ATM/MON/atm_iitmesm_seasonal_hist_clim_mon_jjas_ensmean.nc'
ds_anthrop_clim = xr.open_dataset(data_dir).sel(lat=slice(late,lats)).sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon')

data_dir = '/home/cccr/msingh/phd_data_cmip_account/no_anthrop/ATM/MON/atm_iitmesm_seasonal_hist_clim_mon_jjas_ensmean.nc'
ds_no_anthrop_clim = xr.open_dataset(data_dir).sel(lat=slice(late,lats)).sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon')


for i_y,year_ in enumerate(range(1982,2014)):
    #lats, late, lons, lone = 14, 28, 74, 87
    print(year_)
    year = str(year_)
    ds_imd_ = ds_imd.sel(time=slice(year+'-06-01T00:00:00.000000000',year+'-09-30T00:00:00.000000000')).rf.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone)).mean(dim='lat', skipna=True).mean(dim='lon', skipna=True)
    
    #ds_imd_.mean(dim='time').plot()
    ####################################
    data_dir = '/home/cccr/msingh/phd_data_cmip_account/anthrop/ATM/MON/'+str(year_)
    filenames_a = glob.glob(data_dir+'/atm_iitmesm_seasonal_hist_*_mon_jjas.nc')
    if len(filenames_a)==0:
        continue
    #print(len(filenames))
    ens_list = []
    for file in filenames_a:
        ens_list.append(xr.open_dataset(file).pr.sel(lat=slice(late,lats)).sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon'))
    ds_anthrop_ = xr.concat(ens_list, dim='time')
    ds_anthrop_anom = ds_anthrop_.values
    for i_ens in range(len(filenames_a)):
        ds_anthrop_anom[i_ens*4:(i_ens+1)*4] = ds_anthrop_.values[i_ens*4:(i_ens+1)*4] - ds_anthrop_clim.pr.values
        #print('check')

    ####################################
    
    ####################################
    data_dir = '/home/cccr/msingh/phd_data_cmip_account/no_anthrop/ATM/MON/'+str(year_)
    filenames_na = glob.glob(data_dir+'/atm_iitmesm_seasonal_hist_*_mon_jjas.nc')
    if len(filenames_na)==0:
        continue
    ens_list = []
    #print(len(filenames))
    for file in filenames_na:
        ens_list.append(xr.open_dataset(file).pr.sel(lat=slice(late,lats)).sel(lon=slice(lons,lone)).mean(dim='lat').mean(dim='lon'))
    ds_no_anthrop_ = xr.concat(ens_list, dim='time')

    ds_no_anthrop_anom = ds_no_anthrop_.values
    for i_ens in range(len(filenames_na)):
        ds_no_anthrop_anom[i_ens*4:(i_ens+1)*4] = ds_no_anthrop_.values[i_ens*4:(i_ens+1)*4] - ds_no_anthrop_clim.pr.values

    ####################################
    
    #print(np.sum(np.isnan(ds_imd_.values[0,:,:])), ds_imd_.values[0,:,:].shape[0]* ds_imd_.values[0,:,:].shape[1])
    ds_imd_o = ds_imd_
    ds_imd_jjas_o = ds_imd_o.resample(time='m').mean()
    #print(np.sum(np.isnan(ds_imd_jjas_o.values[0,:,:])), ds_imd_jjas_o.values[0,:,:].shape[0]* ds_imd_jjas_o.values[0,:,:].shape[1])
    #lats, late, lons, lone = 14, 28, 74, 87
    ds_anthrop_mt = ds_anthrop_anom#.sel(lat=slice(late,lats)).sel(lon=slice(lons,lone))
    ds_no_anthrop_mt = ds_no_anthrop_anom#.sel(lat=slice(late,lats)).sel(lon=slice(lons,lone))
    #print(lats,late,lons,lone)
    #print(ds_anthrop_anom.shape, ds_no_anthrop_anom.shape)
    ds_imd_mt = ds_imd_jjas_o#.sel(lat=slice(lats,late)).sel(lon=slice(lons,lone))
    #print(ds_imd_mt)
    #print(np.sum(np.isnan(ds_imd_mt.values[0,:,:])), ds_imd_mt.values[0,:,:].shape[0]* ds_imd_mt.values[0,:,:].shape[1])
    #print(ds_imd_mt.values[0,0,0],  ds_anthrop_mt[0,0,0]) - this is fine
    
    ds_imd_mt_rep_anthrop = np.repeat([ds_imd_mt.values], [len(filenames_a)], axis=0)
    ds_imd_mt_rep_anthrop_ = ds_anthrop_mt.copy()#.values
    #print(ds_imd_mt.values[0],  ds_anthrop_mt[0])
    
    for i_rep in range(len(filenames_a)):
        ds_imd_mt_rep_anthrop_[i_rep*4:(i_rep+1)*4,] = ds_imd_mt_rep_anthrop[i_rep,:]
        
    #print(ds_imd_mt.values[0,0,0],  ds_anthrop_mt[0,0,0])
    ds_imd_mt_rep_no_anthrop = np.repeat([ds_imd_mt.values], [len(filenames_na)], axis=0)
    ds_imd_mt_rep_no_anthrop_ = ds_no_anthrop_mt.copy()#.values
    for i_rep in range(len(filenames_na)):
        ds_imd_mt_rep_no_anthrop_[i_rep*4:(i_rep+1)*4, ] = ds_imd_mt_rep_no_anthrop[i_rep,:]

    #print(ds_imd_mt.values[0,0,0], ds_anthrop_mt[0,0,0])

    #print(y_actual, y_predicted)
    y_actual =  ds_imd_mt_rep_anthrop_.flatten()# ds_imd_mt.values.flatten()
    y_predicted = ds_anthrop_mt.flatten()*86400
    #print(y_actual.shape, y_predicted.shape)
    corr_anthrop =  ma.corrcoef(ma.masked_invalid(y_actual), ma.masked_invalid(y_predicted))

    y_actual = ds_imd_mt_rep_no_anthrop_.flatten()# ds_imd_mt.values.flatten()
    y_predicted = ds_no_anthrop_mt.flatten()*86400
    #print(y_actual.shape, y_predicted.shape)
    corr_no_anthrop = ma.corrcoef(ma.masked_invalid(y_actual), ma.masked_invalid(y_predicted))


    
    corr_anthrop_[i_y] = corr_anthrop[0,1]
    corr_no_anthrop_[i_y] = corr_no_anthrop[0,1]
