def soilm_memory(ds_soilm):
    threshold    = 1./np.exp(1.)
    soilm = ds_soilm.values
    smemory = np.zeros_like((ds_piClim_control_MPIESM_r1i1p1f1_mrso.mrso.values[0,:,:]))
    ntim = 150 # ds_piClim_control_MPIESM_r1i1p1f1_mrso.mrso.values.shape[0]
    correlation_ = np.zeros((ntim, \
                        ds_piClim_control_MPIESM_r1i1p1f1_mrso.mrso.shape[1], \
                        ds_piClim_control_MPIESM_r1i1p1f1_mrso.mrso.shape[2]))    
    for tt in range(2,ntim):
        soilm_lagged  =  soilm[tt:,:,:]
        times = ds_piClim_control_MPIESM_r1i1p1f1_mrso.time.values[tt:]
        lats = ds_piClim_control_MPIESM_r1i1p1f1_mrso.lat.values
        lons = ds_piClim_control_MPIESM_r1i1p1f1_mrso.lon.values
        ds = xr.Dataset({
        'soilm': xr.DataArray(
                    data   = soilm[:-tt],   # enter data here
                    dims   = ['time', 'lat', 'lon'],
                    coords = {'time': times, 'lat':lats, 'lon':lons},
                    ),
         'soilm_lagged': xr.DataArray(
                    data   = soilm_lagged,   # enter data here
                    dims   = ['time', 'lat', 'lon'],
                    coords = {'time': times, 'lat':lats, 'lon':lons},

                    )
                },
        )
        print('lag = ', tt)
        x = ds['soilm']
        y = ds['soilm_lagged']
        n = y.notnull().sum(dim='time')
        xmean = x.mean(axis=0)
        ymean = y.mean(axis=0)
        xstd  = x.std(axis=0)
        ystd  = y.std(axis=0)

        #4. Compute covariance along time axis
        cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

        #5. Compute correlation along time axis
        cor   = cov/(xstd*ystd)
        correlation_[tt-2,:,:] = cor.values
    for i_lat in range(correlation_.shape[1]):
        for j_lon in range(correlation_.shape[2]):
            idx = np.where(correlation_[:,i_lat, j_lon] < 1/np.exp(1.))[0]
            #print(idx)
            #print(len(idx))
            #print(np.sum(np.isnan(soilm[:,i_lat,j_lon])))
            if len(idx)==0:
                smemory[i_lat, j_lon] = np.nan
            elif len(idx)==2:
                smemory[i_lat, j_lon] = np.nan
            else:
                print(idx)
                smemory[i_lat, j_lon] = idx[0]+1
    ds = xr.Dataset({
    'smemory': xr.DataArray(
                data   = smemory,   # enter data here
                dims   = [ 'lat', 'lon'],
                coords = {'lat':lats, 'lon':lons},
                )
            },
    )
            
    return ds
