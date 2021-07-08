# Source: https://unidata.github.io/MetPy/latest/tutorials/xarray_tutorial.html

data_at_point = data.metpy.sel(
    time1='2017-09-05 12:00',
    latitude=40 * units.degrees_north,
    longitude=260 * units.degrees_east
)
dewpoint = mpcalc.dewpoint_from_relative_humidity(
    data_at_point['Temperature_isobaric'],
    data_at_point['Relative_humidity_isobaric']
)
cape, cin = mpcalc.surface_based_cape_cin(
    data_at_point['isobaric3'],
    data_at_point['Temperature_isobaric'],
    dewpoint
)
cape

