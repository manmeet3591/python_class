import os
import pandas as pd
#os.system(" acccmip6 -o S -m CanESM5 -e historical -v pr -f day"
models = pd.read_csv("/lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/CMIP6/models_esgf_cmip_historical_pr.txt", header=None).values
print(models.shape)
for i in range(models.shape[0]):
    print(models[i,0])
    os.system("/lus/dal/cccr_rnd/manmeet/anaconda3/envs/py36/bin/acccmip6 -o S -m "+models[i,0]+" -e historical -v pr -f day")
    os.system("/lus/dal/cccr_rnd/manmeet/anaconda3/envs/py36/bin/acccmip6 -o D -m "+models[i,0]+" -e historical -v pr -f day")
