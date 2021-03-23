# References
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011GL048268
# http://www.coupling-metrics.com/terrestrial-coupling-index/

# Calculate the standard deviation of soil moisture.
# Calculate the slope between soil moisture and surface flux (either sensible or latent) with soil moisture as the independent variable.
# Multiple the standard deviation by the slope returning the terrestrial coupling index with the same units as the surface flux variable.

def linear_trend(x, y):
    if np.sum(np.isnan(x))>0:
        pf = np.empty((x.shape[0]))
    else:
    #idx = np.isfinite(x) & np.isfinite(y)
        pf = np.polyfit(x, y, 1)
    return xr.DataArray(pf[0])

def compute_tci(sm, lh):
    sm_std = sm.std(dim='time').plot()
    slopes = xr.apply_ufunc(linear_trend,
                            sm, lh,
                            vectorize=True,
                            input_core_dims=[['time'], ['time']],# reduce along 'plev'
                            )
    tci = slopes*sm_std
#mrsos_control[0].std(dim='time').plot()
