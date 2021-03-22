#!/bin/bash
#
#for year in {1990..2008}
#do
#python download_era5.py $year evaporation_from_vegetation_transpiration
#done
#
#for year in {2010..2020}
#do
#python download_era5.py $year evaporation_from_vegetation_transpiration
#done
#
for year in {1981..2020}
do
python download_era5.py $year total_evaporation
python download_era5.py $year surface_latent_heat_flux
python download_era5.py $year surface_sensible_heat_flux
python download_era5.py $year volumetric_soil_water_layer_1
python download_era5.py $year forecast_albedo
done
