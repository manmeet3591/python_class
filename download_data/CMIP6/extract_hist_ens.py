#!/usr/bin/env python
# coding: utf-8

# In[1]:


from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import zarr
import gcsfs
import seaborn as sns
import os 

xr.set_options(display_style='html')
#get_ipython().run_line_magic('matplotlib', 'inline')
#get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")


# In[2]:


df = pd.read_csv('https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv')
print(df.head())
#
#
## In[3]:
#
#
#df_ta = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'tas' & experiment_id == 'historical'")
#print(df_ta)
#
#
## In[4]:
#
#
#df_pr = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'pr' & experiment_id == 'historical'")
#print(df_pr)
#
#
## In[5]:
#
#
#df_tasmax = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'tasmax' & experiment_id == 'historical'")
#print(df_tasmax)
#
#
## In[6]:
#
#
#df_tasmin = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'tasmin' & experiment_id == 'historical'")
#print(df_tasmin)
#

df_uas = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'uas' & experiment_id == 'historical'")
print(df_uas)
df_vas = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'vas' & experiment_id == 'historical'")
print(df_vas)
df_hurs = df.query("activity_id=='CMIP' & table_id == 'day' & variable_id == 'hurs' & experiment_id == 'historical'")
print(df_vas)
# In[7]:


#df_pr['institution_id'].value_counts().to_csv('pr_historical_data_available_unique_models.csv')
#df_ta['institution_id'].value_counts().to_csv('tas_historical_data_available_unique_models.csv')
#df_tasmax['institution_id'].value_counts().to_csv('tasmax_historical_data_available_unique_models.csv')
#df_tasmin['institution_id'].value_counts().to_csv('tasmin_historical_data_available_unique_models.csv')
df_uas['institution_id'].value_counts().to_csv('uas_historical_data_available_unique_models.csv')
df_vas['institution_id'].value_counts().to_csv('vas_historical_data_available_unique_models.csv')
df_hurs['institution_id'].value_counts().to_csv('hurs_historical_data_available_unique_models.csv')


# In[20]:


var_ = ['pr', 'tas', 'tasmax', 'tasmin', 'uas', 'vas', 'hurs']
#var_ = ['uas', 'vas']
#var_ = ['hurs']
for var in var_:
    #print(var, type(var))
    unique_models = pd.read_csv(var+'_historical_data_available_unique_models.csv')
    #unique_models.iloc[:,0]
    for i_model,model in enumerate(unique_models.iloc[:,0]):
        #tas_model_ = []
        #print('test')
        #print(type(model))
        #print("table_id == 'day' & variable_id == 'tas' & experiment_id == 'historical' & institution_id =='"+model+"'")
        query_ = "activity_id=='CMIP' & table_id == 'day' & variable_id == '"+var+"' & experiment_id == 'historical' & institution_id =='"+model+"'"
        df_ta = df.query(query_)
        #print('test')
        filename_ = var+'_historical_'+str(model)+".nc"
#        filename_ = filename_[:-3]+'_ens_'+str(ens+1)+'.nc'
        for ens in range(unique_models.iloc[i_model,1]):
            filename_ = var+'_historical_'+str(model)+'_ens_'+str(ens+1)+'.nc'
            count = []
            #print(unique_models.iloc[i_model,1])
            if not os.path.isfile(filename_):
            #for ens in range(1):
            # this only needs to be created once
                
                gcs = gcsfs.GCSFileSystem(token='anon')
                #print('test')
                # get the path to a specific zarr store (the first one from the dataframe above)
                zstore = df_ta.zstore.values[ens] # The index here controls everything, the number

                # create a mutable-mapping-style interface to the store
                mapper = gcs.get_mapper(zstore)

                # open it using xarray and zarr
                ds_ta = xr.open_zarr(mapper, consolidated=True)
                #print(ds_ta.sel(lat=slice(29,32)).sel(lon=slice(261,264)))
                
                try:
                    #tas_hist_ = ds_ta[var].sel(lat=30.3134, lon= 360-97.7666, method='nearest').sel(time=slice('2021','2098'))
                    #print('test 1') #tas_hist_ = ds_ta[var].sel(lat=30.3134, lon= 360-97.7666, method='nearest').sel(time=slice('2021','2098'))
                    tas_hist_ = ds_ta[var].sel(lat=slice(29,32)).sel(lon=slice(261,264)).sel(time=slice('1950','2014'))
                    #print(tas_hist_.values.shape[0]>0) 

                    if tas_hist_.values.shape[0]>0:
                        count.append(ens)
                    #print('After ', tas_hist_.values.shape) 

                except:
                    continue
                
                #print(tas_hist_.time.values.shape[0])
                #print(type(tas_hist_.time.values[0]).__name__)
                #no_of_days=5475
                #print(count)
                
                print("Writing ", filename_)
                out_data = tas_hist_
                out_data[var+'_'] = (('time', 'lat', 'lon'), np.asarray(tas_hist_))
                out_data[var+'_'].to_netcdf(filename_)
#                if tas_hist_.values.shape[0]>0:
#                    if ens==count[0]:
#                        if type(tas_hist_.time.values[0]).__name__=='DatetimeNoLeap':
#                            no_of_days=28470
#                        elif type(tas_hist_.time.values[0]).__name__=='Datetime360Day':
#                            no_of_days=28080
#                        else: 
#                            no_of_days=28489
#                        out_data = tas_hist_
#                else:
#                    continue
#                print('test', no_of_days, tas_hist_.shape)
#                if tas_hist_.values.shape[0]==no_of_days:
#		    # There is a difference of resolution within the same model for some cases
#                    if ens==0 and tas_hist_.shape[1]*tas_hist_.shape[2]<5: 
#                        tas_model_.append(tas_hist_.values)
#			
#                    if ens>0: 
#                        if tas_hist_.shape[1]*tas_hist_.shape[2]<5: 
##                            if tas_model_[0].shape==tas_hist_.shape: 
##                            tas_model_.append(tas_hist_.values) if tas_model_[0].shape==tas_hist_.shape else continue
#                            if tas_model_[0].shape==tas_hist_.shape: tas_model_.append(tas_hist_.values)
#                    if tas_hist_.shape[1]*tas_hist_.shape[2]>4 and not os.path.isfile(filename_[:-3]+'_ens_'+str(ens+1)+'.nc'): 
#                        out_data[var+'_'] = (('time', 'lat', 'lon'), np.asarray(tas_hist_))
#                        out_data[var+'_'].to_netcdf(filename_[:-3]+'_ens_'+str(ens+1)+'.nc')
#                        print('Writing '+filename_[:-3]+'_ens_'+str(ens+1)+'.nc')
#                print('ensemble number: ',ens)
#                ##############################################################################
#            #print('test ', len(tas_model_))
#            if len(tas_model_)>0:
#
##                print(tas_model_[0].shape[0], tas_model_[0].shape[1],tas_model_[0].shape[2])
#                #tas_model__ = np.zeros((tas_model_[0].shape[0], tas_model_[0].shape[1],tas_model_[0].shape[2]))
#                #print(np.nanmean(np.asarray(tas_model_),axis=0).shape)
#                #tas_model__ = np.nanmean(np.asarray(tas_model_),axis=0)
#                #print(tas_model__.shape)
#
#                if tas_hist_.shape[1]*tas_hist_.shape[2]<5:
#                    out_data[var+'_'] = (('time', 'lat', 'lon'), np.nanmean(np.asarray(tas_model_),axis=0))
#                    out_data[var+'_'].to_netcdf(filename_)
#                else:
#                    continue
#                        #out_data[var+'_'] = (('time', 'lat', 'lon'), np.asarray(tas_model_[i_ens]))
#                        #out_data[var+'_'].to_netcdf(filename_[:-3]+'_ens_'+str(i_ens+1)+'.nc')
#                



    #print(df_ta)
# unique_models.iloc[:,0]


# In[ ]:

#
#count[0]
#
#
## In[ ]:
#
#
#tas_model_.shape
#

# In[ ]:





# In[ ]:





# In[ ]:




