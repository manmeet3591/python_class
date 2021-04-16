# TCI for control
# mrsos_control[0].shape
# mrsos_std = np.zeros((mrsos_control[0].shape[1], mrsos_control[0].shape[2]))
# for i_lat in range(mrsos_control[0].shape[1]):
#     for j_lon in range(mrsos_control)
# https://stackoverflow.com/questions/58546999/calculate-correlation-in-xarray-with-missing-data
def linear_trend(x, y):
    #print(x.shape)
#     if np.sum(np.isnan(x))>0:
#         pf = np.empty((x.shape[0]))
#     else:
#         try:
    idx = np.isfinite(x) & np.isfinite(y)
    #print(x.shape)
    pf = np.polyfit(x[idx], y[idx], 1)
#         except:
#             pf = np.empty((x.shape[0]))
    return xr.DataArray(pf[0])

def compute_tci(sm, lh):
    print(sm.shape, lh.shape)
    sm_std = sm.std(dim='time')
    slopes = np.zeros_like(sm_std.values)
#     slopes = np.zeros_like(sm_std)
#     for i_tci in range(sm.shape[1]):
#         for j_tci in range(sm.shape[2]):
#             mask = ~np.isnan(sm) & ~np.isnan(lh)
#             slope, intercept, r_value, p_value, std_err = stats.linregress(varx[mask], vary[mask])
#             slopes[i_tci,j_tci] = slope
#     for i_tci in range(sm.values.shape[1]):
#         for j_tci in range(sm.values.shape[2]):
#             mask = ~np.isnan(sm.values[:,i_tci,j_tci]) & ~np.isnan(lh.values[:,i_tci,j_tci])
#             if np.sum(np.isnan(sm.values[:,i_tci,j_tci]))>=sm.values.shape[0]-2:
#                 slope=np.nan
#             else:
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(sm.values[:,i_tci,j_tci][mask], lh.values[:,i_tci,j_tci][mask])
#             slopes[i_tci,j_tci] = slope
            #print(i_tci, j_tci)
#     slopes = xr.apply_ufunc(linear_trend,
#                             sm, lh,
#                             vectorize=True,
#                             input_core_dims=[['time'], ['time']],# reduce along 'plev'
#                             )
    x = sm
    y = lh
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)
    slopes     = cov/(xstd**2)
    intercept = ymean - xmean*slope
    sm_std['slopes'] = (('lat', 'lon'), slopes)
    tci = sm_std.slopes*sm_std
    return tci
#mrsos_control[0].std(dim='time').plot()

mrsos_control = ds_piClim_control_MPIESM_r1i1p1f1_mrso.mrso
hfls_control  = ds_piClim_control_MPIESM_r1i1p1f1_hfls.hfls
tci_control = compute_tci(mrsos_control.sel(time=mrsos_control.time.dt.month.isin([5,6, 7, 8,9])), \
                                  hfls_aer.sel(time=mrsos_control.time.dt.month.isin([5,6, 7, 8,9])))
