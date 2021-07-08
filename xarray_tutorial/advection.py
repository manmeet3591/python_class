import metpy.calc as mpcalc
temp_850 = temp[lev_850]
u_wind_850 = u_wind[lev_850]
v_wind_850 = v_wind[lev_850]

# Use helper function defined above to calculate distance
# between lat/lon grid points
dx, dy = mpcalc.lat_lon_grid_deltas(lon_var, lat_var)

# Calculate temperature advection using metpy function
adv = mpcalc.advection(temp_850 * units.kelvin, [u_wind_850, v_wind_850],
                       (dx, dy), dim_order='yx') * units('K/sec')
