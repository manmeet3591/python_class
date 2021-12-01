ds = u.to_dataset(name = 'uas')
ds['vas'] = v

ds.plot.quiver('lon', 'lat', 'uas', 'vas', ax=ax[0,0])
