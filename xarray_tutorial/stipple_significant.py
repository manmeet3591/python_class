fig = plt.figure(figsize=(11,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8],projection=ccrs.PlateCarree())
a =  ds_merra2_mam.dust_trend_mam.sel(lon=slice(60,100), lat=slice(5,40)).plot(ax=ax, \
                                                                          vmin=-0.01,\
                                                                          vmax=0.01,extend='both', \
                                                                          cmap='twilight_shifted')
ds_merra2_mam['dust_p_value_manm_'] = ds_merra2_mam['dust_p_value_mam']<0.05
#ds_merra2_mam['dust_p_value_manm_'].sel(lon=slice(60,100), lat=slice(5,40)).plot.scatter(ax=ax)

[m,n] = np.where(ds_merra2_mam['dust_p_value_mam'].sel(lon=slice(60,100), lat=slice(5,40)).values<0.05 )
z1 = np.zeros(ds_merra2_mam['dust_p_value_mam'].sel(lon=slice(60,100), lat=slice(5,40)).values.shape)
z1[m, n] = 99
x = ds_merra2_mam.sel(lon=slice(60,100), lat=slice(5,40)).lon.values
y = ds_merra2_mam.sel(lon=slice(60,100), lat=slice(5,40)).lat.values

cs = ax.contourf(x, y, z1 ,1 , hatches=['', '..'],  alpha=0)

ax.coastlines(resolution='10m')
plt.title('Dust AOD(550nm) MAM trend from MERRA2 reanalysis 2002-2020')
ax_ = [ax]
for i in range(1):
    gl = ax_[i].gridlines()
    gl.xlabels_bottom = True
    gl.ylabels_left = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 'medium'}
    gl.ylabel_style = {'size': 'medium'}
#plt.savefig('merra2_mam_dust_22_34n_65_90e_2002_2014_spatial_reproduce_vinoj.png')
