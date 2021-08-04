# python_class
This repository contains all the scripts used for the python class for JRFs at IITM

Some of the tutorials to follow
1. http://pure.iiasa.ac.at/id/eprint/14952/1/xarray-tutorial-egu2017-answers.pdf
2. https://rabernat.github.io/research_computing/xarray.html
3. Ocean Data Analysis https://currents.soest.hawaii.edu/ocn_data_analysis/exercise_data.html#id1
4. Parallelization http://xarray.pydata.org/en/stable/dask.html
5. Satellite Data Analyis https://github.com/nansencenter/nansat-lectures
6. https://github.com/NCAR/CESM_postprocessing CESM Postprocessing
7. https://github.com/NCAR/PyCect This repo is used to compare the results of a set of new CAM simulations against the accepted ensemble
8. https://github.com/nichannah/ocean-regrid Regrid ocean reanalysis data from normal to tripolar grids
9. https://github.com/jswhit/gfstonc Read GFS sigma and sfc files in python
10. f2py
11. Pandas
12. https://github.com/tmiyachi/data2gfs Make python version of this using f2py
13. Shallow water equation model using pyspharm https://github.com/jswhit/pyspharm and https://www.aosc.umd.edu/~dkleist/docs/shtns/doc/html/shallow_water_8py-example.html explaining the code
14. Scientific Computing Lectures https://github.com/jrjohansson/scientific-python-lectures
15. Geopandas satellite data analysis https://towardsdatascience.com/satellite-imagery-access-and-analysis-in-python-jupyter-notebooks-387971ece84b
16. Rasterio https://medium.com/analytics-vidhya/satellite-imagery-analysis-with-python-3f8ccf8a7c32
17. Eo-learn https://medium.com/dataseries/satellite-imagery-analysis-with-python-ii-8001e5c41a52
18. Satpy
19. Use of Landsat and Sentinel datasets
20. Pyunicorn
21. Keras, tensorflow, pytorch, django, theano, scikit-learn, theano, bokeh, pandas, seaborn, bokeh, plotly, scrapy, 
22. Python tutorial https://carpentrieslab.github.io/python-aos-lesson/ plotting CMIP data - highlight
23. Python for oceanography http://www.soest.hawaii.edu/oceanography/courses/OCN681/python.html
24. Python tools for oceanography https://pyoceans.github.io/sea-py/
25. Python Land Surface Modelling https://www.geosci-model-dev.net/12/2781/2019/
26. Python hydrology tools https://github.com/raoulcollenteur/Python-Hydrology-Tools
27. Docker
28. Python and GIS https://automating-gis-processes.github.io/CSC18/lessons/L1/overview.html
29. https://automating-gis-processes.github.io/2016/
30. https://geohackweek.github.io/raster/
31. https://github.com/pangeo-data/pangeo 
32. https://github.com/pangeo-data/awesome-open-climate-science
33. https://uwescience.github.io/sat-image-analysis/resources.html
34. Radar data analysis https://data.world/datasets/radar https://arm-doe.github.io/pyart/ https://docs.wradlib.org/
35. https://www.earthdatascience.org/courses/use-data-open-source-python/multispectral-remote-sensing/landsat-in-Python/
36. Deep Learning on Satellite Imagery https://github.com/robmarkcole/satellite-image-deep-learning
37. Google Earth Engine https://sites.google.com/view/eeindia-advanced-summit/summit-resources
38. https://geohackweek.github.io/GEE-Python-API/
39. https://github.com/google/earthengine-api/tree/master/python/examples/ipynb
40. http://www.jerico-ri.eu/download/summer%20school%20-%20the%20netherlands/Genna%20Donchyts%20-%20GEE%20Training.pdf
41. https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
42. Installing Google Earth Engine and requesting access https://github.com/google/earthengine-api/issues/27
43. https://github.com/giswqs/earthengine-py-notebooks
44. Google Earth Engine image to numpy https://mygeoblog.com/2019/08/21/google-earth-engine-to-numpy/
45. Stippling to show statistical significance https://github.com/bradyrx/esmtools/issues/13
46. Resampling from swath to grid https://github.com/TerraFusion/pytaf
47. Making a docker container for data science https://towardsdatascience.com/docker-for-data-scientists-5732501f0ba4
48. Docker commands:

Run interactively: docker run -it  manmeet3591/dl:iitm:latest

Install the necessary libraries

Open a new terminal and do docker images to see the id and run the following command

$ docker tag id_ manmeet3591/dl_iitm:v2

$ docker push manmeet3591/dl_iitm:v2

Projects for the class

https://docs.google.com/spreadsheets/d/1m2ZIJ_To8IbE18Teb70a7BVZg0o29sOM6rlgFkE2b3E/edit#gid=0

https://docs.google.com/document/d/12h9bcIdBPJUFc_fJssJe8hVBzledq2Dtk5-9OpKHbfg/edit 

49. Homogenous regions India shape files: https://github.com/Cassimsannan/Shapefiles

50. Download CMIP6 data: https://github.com/TaufiqHassan/acccmip6

51. Download MSWEP data from Google drive: 

Setup rclone: https://www.youtube.com/watch?v=vPs9K_VC-lg

52. Run jupyter notebook from docker container 

docker run --rm -it --entrypoint bash -p 8891:8891 manmeet3591/tensortrade

Inside the container 
jupyter-notebook --ip 0.0.0.0 --port=8891 --no-browser --allow-root &

In the browser 
http://localhost:8891/

$ rclone sync -v --exclude 3hourly/ --drive-shared-with-me GoogleDrive:/MSWEP_V280 /lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/mswep/.

53. Create any number of subplots matplotlib

$ fig,ax = plt.subplots(ncols=2,nrows=4, figsize=(11.69,8.27), subplot_kw={'projection': ccrs.PlateCarree()})

54. Google Earth Engine timelapse gif generator: https://9611d0317f71.ngrok.io/voila/render/timelapse.ipynb

55. Handling expver dimension in a netcdf file downloaded as ERA5 data 

 ds.reduce(np.nansum, 'expver') Solution from marco venturini https://confluence.ecmwf.int/pages/viewpage.action?pageId=173385064
 
 56. GeoTIFF to netcdf and exporting data from Google Earth Engine https://medium.com/@wenzhao.li1989/nco-translate-geotiff-files-exported-from-gee-to-a-netcdf-file-with-correct-time-dimension-ce97a8f3043f

57. Make pipeline to avoid test data leaking into train https://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html

58. Potential evapotranspiration (PET) from netcdf file https://climate-indices.readthedocs.io/en/latest/#

59. t-distributed Stochastic Neighbourhood Embedding (tSNE) versus PCA https://stats.stackexchange.com/questions/238538/are-there-cases-where-pca-is-more-suitable-than-t-sne

60. Stationarity of time series: https://towardsdatascience.com/stationarity-in-time-series-analysis-90c94f27322

61. SARIMAX model: https://towardsdatascience.com/end-to-end-time-series-analysis-and-forecasting-a-trio-of-sarimax-lstm-and-prophet-part-1-306367e57db8

62. Prevent kaggle from disconnecting https://stackoverflow.com/questions/57113226/how-to-prevent-google-colab-from-disconnecting

function ClickConnect(){
    console.log("Working"); 
    document.querySelector("colab-toolbar-button#connect").click() 
}
setInterval(ClickConnect,60000)

63. Solving NVIDIA driver installation issues https://stackoverflow.com/questions/42984743/nvidia-smi-has-failed-because-it-couldnt-communicate-with-the-nvidia-driver/51113428#51113428

64. Installing NVIDIA drivers https://www.itzgeek.com/post/how-to-install-nvidia-drivers-on-ubuntu-20-04-ubuntu-18-04.html

65. Install cuda https://www.tensorflow.org/install/gpu

66. bashrc commands for cuda

export CUDA_HOME=/usr/local/cuda-11.0
export LD_LIBRARY_PATH=/usr/local/cuda-11.0/lib64:/usr/local/cuda-11.0/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/cuda-11.0/bin:$PATH



Troubleshooting

1. Continue in outer loop using multi-loops https://stackoverflow.com/questions/14829640/how-to-continue-in-nested-loops-in-python

2. Numbering the subplots https://matplotlib.org/3.1.1/gallery/axes_grid1/simple_anchored_artists.html

3. Fortran compilation may sometimes be solved by running the command ulimit -s unlimited

4. There are visualization problems in cartopy if the lon is from 0 to 360 and not from -180 to 180

5. Run docker as a non-root user https://docs.docker.com/engine/install/linux-postinstall/

6. In the first instance of an image sometimes docker hub may deny you to push the image https://stackoverflow.com/questions/41984399/denied-requested-access-to-the-resource-is-denied-docker

7. Numpy to xarray : foo = xr.DataArray(data, coords=[times, locs], dims=["time", "space"])

  data = ds_merra2_jjas.DUSCATAU.sel(time='2002').values[0,:,:]
  lats_ = ds_merra2_jjas.DUSCATAU.sel(time='2002').lat.values
  lons_ = ds_merra2_jjas.DUSCATAU.sel(time='2002').lon.values
  ds_merra2_jjas_new = xr.DataArray(data, coords=[lats_, lons_], dims=["lat", "lon"])

8. Using matplotlib to make map plots
  plt.contourf(ds_merra2_jjas.DUSCATAU.sel(time='2002').lon.values,\
                 ds_merra2_jjas.DUSCATAU.sel(time='2002').lat.values ,\
                ds_merra2_jjas.DUSCATAU.sel(time='2002').values[0,:,:],\
                cmap='bwr')
plt.colorbar()

9. Sometimes xarray plot might show blank, the way to resolve that is select the area and that should work.

10. Pattern correlation formula: https://www.mdpi.com/2073-4441/10/1/28 may use weights as well for the pattern correlation 

For the weights, the following can be followed: https://stackoverflow.com/questions/58881607/calculating-the-cosine-of-latitude-as-weights-for-gridded-data

11. When installing packages otherwise difficult to install like ESMF we can set the compiler environment variables such as CC and FC to force conda to install using that particular compiler. This saves a lot of time and effort. https://stackoverflow.com/questions/59284298/conda-install-c-anaconda-gcc-linux-64-not-being-used  Many build tools such as make and CMake search by default for a compiler named simply gcc, so we set environment variables to point these tools to the correct compiler.

12. When using the isin function with sel we can at present use it only once in a call. Need to instantiate a new variable for doing it twice.

13. Installing PyRQA (Runs only with python 2.7) https://github.com/szhan/pyrqa

    conda install https://anaconda.org/conda-forge/pytools/2017.2/download/linux-64/pytools-2017.2-py27_0.tar.bz2
    
    conda install https://anaconda.org/conda-forge/pyopencl/2018.1.1/download/linux-64/pyopencl-2018.1.1-py27_1.tar.bz2
    
    conda install -c conda-forge pocl
    
    pip install Mako
    
    pip install PyRQA
    
    Even after all this, unable to run pyrqa smoothly. However, this activity ensured that the environment to run pyrqa was perfect. So then clone the github repository and inside the main github repository pyrqa, there is a folder pyrqa. Copy that to your desired location, rename it lets say PYRQA. And use the library as PYRQA.

14. Logging to a remote server without password https://www.thegeekstuff.com/2008/11/3-steps-to-perform-ssh-login-without-password-using-ssh-keygen-ssh-copy-id/





