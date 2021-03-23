# References
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011GL048268
# http://www.coupling-metrics.com/terrestrial-coupling-index/

# Calculate the standard deviation of soil moisture.
# Calculate the slope between soil moisture and surface flux (either sensible or latent) with soil moisture as the independent variable.
# Multiple the standard deviation by the slope returning the terrestrial coupling index with the same units as the surface flux variable.

def compute_tci(sm, lh):
    print(sm.shape, lh.shape)
    sm_std = sm.std(dim='time')
    slopes = np.zeros_like(sm_std.values)

    for i_tci in range(sm.values.shape[1]):
        for j_tci in range(sm.values.shape[2]):
            mask = ~np.isnan(sm.values[:,i_tci,j_tci]) & ~np.isnan(lh.values[:,i_tci,j_tci])
            if np.sum(np.isnan(sm.values[:,i_tci,j_tci]))>=sm.values.shape[0]-2:
                slope=np.nan
            else:
                slope, intercept, r_value, p_value, std_err = stats.linregress(sm.values[:,i_tci,j_tci][mask], lh.values[:,i_tci,j_tci][mask])
            slopes[i_tci,j_tci] = slope
        

    sm_std['slopes'] = (('lat', 'lon'), slopes)
    tci = sm_std.slopes*sm_std
    return tci
