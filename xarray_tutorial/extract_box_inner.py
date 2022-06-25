ds_o = ds.where((ds["lon"] < 60) | (ds["lon"] > 80), drop=True).where((ds["lat"] < 10) | (ds["lat"] > 20), drop=True)
