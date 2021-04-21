import pandas as pd
import numpy as np
import pyet as pyet
def calc_pet_oudin(tas_control):
    pet_control = np.zeros_like((tas_control.values))
    time_pd = pd.to_datetime(tas_control.time.values)
    for i_lat in range(tas_control.shape[1]):
        print(i_lat)
        for j_lon in range(tas_control.shape[2]):

            df = pd.DataFrame({'YYYYMMDD':time_pd, 'data':tas_control.values[:,i_lat,j_lon]-273.15})
            df.set_index('YYYYMMDD', inplace=True)
            et_oudin = pyet.oudin(df['data'], lat=np.radians(tas_control.lat.values[i_lat])).values
            pet_control[:,i_lat,j_lon]= et_oudin
    return pet_control

months = [6,7,8,9]
tas_control    = ds_piClim_control_MPIESM_r1i1p1f1_tas.tas.sel(time=ds_piClim_aer_MPIESM_r1i1p1f1_hfss.time.dt.month.isin(months))
pet_control = np.zeros_like((tas_control.values))
time_pd = pd.to_datetime(tas_control.time.values)
for i_lat in range(tas_control.shape[1]):
    print(i_lat)
    for j_lon in range(tas_control.shape[2]):
        
        df = pd.DataFrame({'YYYYMMDD':time_pd, 'data':tas_control.values[:,i_lat,j_lon]-273.15})
        df.set_index('YYYYMMDD', inplace=True)
        et_oudin = pyet.oudin(df['data'], lat=np.radians(tas_control.lat.values[i_lat])).values
        pet_control[:,i_lat,j_lon]= et_oudin
