# https://wiki.mpimet.mpg.de/doku.php?id=analysis:pot_pourri:processing_binary_data
import gzip
import xarray as xr
import numpy as np
import pandas as pd
import os

#filename1='/home/cccr/msingh/satellite_data/hokusai.eorc.jaxa.jp/standard/v6/hourly/2002/01/01/gsmap_rnl.20020101.0000.v6.5133.0.dat.gz'


# In[3]:


#gz = gzip.GzipFile(filename1,'rb')
#dd=np.frombuffer(gz.read(),dtype=np.float32)
#pre=dd.reshape((1200,3600))


#lon=np.linspace(0.05,359.95,3600)
#lat=np.linspace(59.95,-59.95,1200)



#unit_conversion = 1./3600.
#gsm_out = '/home/cccr/msingh/satellite_data/GSMaP.v7_0.10deg_2001.nc'
gsm_out = '/home/cccr/msingh/satellite_data/hokusai.eorc.jaxa.jp/standard/v6/hourly/test.nc'
year = '2002'
if (os.path.isfile(gsm_out)): 
	print (gsm_out+' exists, not overwriting')
else:
#	time= pd.date_range(ym+'2001-01-01','2001-12-31', freq='1H')
	#time= pd.date_range('2001-01-01','2001-12-31 23:00:00', freq='1H')
	time= pd.date_range('2002-02-01','2002-02-01 23:00:00', freq='1H')
	lat    = xr.Coordinate('lat', np.linspace(59.95,-59.95,1200) , attrs={'long_name':'latitude','standard_name':'latitude','units': "degrees_north"})
	lon    = xr.Coordinate('lon', np.linspace(0.05,359.95,3600) , attrs={'long_name':'longitude','standard_name':'longitude','units': "degrees_east"})
	da     = np.ndarray(shape=(time.size,lat.size,lon.size))
	
	path = '/home/cccr/msingh/satellite_data/hokusai.eorc.jaxa.jp/standard/v6/hourly/2002/'
	for i in np.arange(time.size):
		gsm_in = path + time[i].strftime("%m/%d")+ '/gsmap_rnl.2002'+time[i].strftime("%m%d.%H")+'00.v6.5133.0.dat.gz'#  gsmap_mvk.2016'+time[i].strftime("%m%d.%H")+'00.v7.0001.0.dat'
		if (os.path.isfile(gsm_in)):
			gz = gzip.GzipFile(gsm_in,'rb')
			dd = np.frombuffer(gz.read(),dtype=np.float32)
			data = dd.reshape((1200,3600))
			#with open(gsm_in,'rb') as f:
		#		data = np.fromfile(f, dtype=np.float32, count = lon.size*lat.size)                
			#d2   = pd.DataFrame(np.reshape(data,(lat.size,lon.size)), lat, lon)
			d2   = pd.DataFrame(data, lat, lon)
			da[i,:,:] = d2.where(d2 != -99.0)# * unit_conversion
			print("read successfull ", d2, d2.shape)
		else:
			print(gsm_in, " file not found: filing nan")
			da[i,:,:] = np.nan
	
	da[da<0] = np.nan 
	ds     = xr.Dataset({'pr': (['time', 'lat', 'lon'],  da.astype('float32'))}, coords={'time': time,'lat': lat,'lon': lon})
	
	ds.pr.attrs['long_name']    ='precipitation'
	ds.pr.attrs['standard_name']='precipitation_flux'
	ds.pr.attrs['units']        ='kg m-2 s-1'
	ds.to_netcdf(gsm_out)  
