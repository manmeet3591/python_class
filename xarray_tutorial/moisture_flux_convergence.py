# https://earthscience.stackexchange.com/questions/9634/how-to-compute-du-dx-and-dv-dy-in-moisture-flux-convergence  

import metpy.calc as mpcalc
q_850 = q[lev_850]
u_wind_850 = u_wind[lev_850]
v_wind_850 = v_wind[lev_850]

# Use helper function defined above to calculate distance
# between lat/lon grid points
dx, dy = mpcalc.lat_lon_grid_deltas(lon_var, lat_var)

# Calculate temperature advection using metpy function
adv = mpcalc.advection(q_850 * units.kelvin, [u_wind_850, v_wind_850],
                       (dx, dy), dim_order='yx') * units('K/sec')
