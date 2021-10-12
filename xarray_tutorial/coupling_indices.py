def compute_distribution(v):
    """
    v: vector de valores enteros
    devuelve un diccionario con la probabilidad de cada valor
    computado como la frecuencia de ocurrencia
    """
    d= defaultdict(int)
    for e in v: d[e]+=1
    s= float(sum(d.values()))
    return dict((k, v/s) for k, v in d.items())
  
def entropy(y):
    """
    Computa la entropia de un vector discreto
    """
    # P(Y)
    Py= compute_distribution(y)
    res=0.0
    for k, v in Py.items():
        res+=v*log2(v)
    return -res
  
  def conditional_entropy(x, y):
    """
    x: vector de numeros reales
    y: vector de numeros enteros
    devuelve H(Y|X)
    """
    # discretizacion de X 
    #print(int(x.size/10))
    # https://stats.stackexchange.com/questions/179674/number-of-bins-when-computing-mutual-information#:~:text=3%20Answers&text=There%20is%20no%20best%20number,on%20histograms%20have%20been%20proposed.
    # Choosing number of bins
    #hx, bx= histogram(x, bins=int(x.size/10), density=True)
    hx, bx= histogram(x, bins=int(np.sqrt(x.size/5)),density=True)

    Py= compute_distribution(y)
    Px= compute_distribution(digitize(x,bx))

    res= 0
    for ey in set(y):
        # P(X | Y)
        x1= x[y==ey]
        condPxy= compute_distribution(digitize(x1,bx))

        for k, v in condPxy.items():
            res+= (v*Py[ey]*(log2(Px[k]) - log2(v*Py[ey])))
            
def mutual_information(x,y):
    return entropy(y) - conditional_entropy(x,y)

#https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf
def normalized_mutual_information(x,y):
    return (2* mutual_information(x,y))/(entropy(x)+entropy(y))
  
from copent import transent
from pandas import read_csv
import numpy as np

url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00381/PRSA_data_2010.1.1-2014.12.31.csv"
prsa2010 = read_csv(url)
# index: 5(PM2.5),6(Dew Point),7(Temperature),8(Pressure),10(Cumulative Wind Speed)
data = prsa2010.iloc[2200:2700,[5,8]].values

te = np.zeros(24)
for lag in range(1,2):
    te[lag-1] = transent(data[:,0],data[:,1],1)
    print(data[:,0].shape)
    str_ = "TE from pressure to PM2.5 at %d hours lag : %f" %(lag,te[lag-1])
    print(str_)

# TE from y to x  
def transfer_entropy(x,y):
    return transent(x,y,40)
  
def multivariate_mutual_information(xs1,xs2,xtar):
    # https://stackoverflow.com/questions/20332750/python-joint-distribution-of-n-variables
    numBins =int(np.sqrt(xs1.shape[0]/5))  # number of bins in each dimension
    #print(xtar)
    xs1 = xs1[np.isfinite(xs1)]
    xs2 = xs2[np.isfinite(xs1)]
    xtar = xtar[np.isfinite(xs1)]
    
    xs2 = xs2[np.isfinite(xs2)]
    xs2 = xs2[np.isfinite(xtar)]
    
    xtar = xtar[np.isfinite(xs2)]
    xtar = xtar[np.isfinite(xtar)]
    #print(xtar)
    data = np.stack((xs1, xs2, xtar)).T
    #data = np.random.randn(100000, 3)  # generate 100000 3-d random data points
    jointProbs3, edges = np.histogramdd(data, bins=numBins)
    jointProbs3 /= jointProbs3.sum()
    jointProbs2, edges = np.histogramdd(data[:,:2], bins=numBins)
    jointProbs2 /= jointProbs2.sum()
    jointProbs1, edges = np.histogramdd(data[:,2:], bins=numBins)
    jointProbs1 /= jointProbs1.sum()
    jointProbs3_ = jointProbs3.copy()
    for i in range(jointProbs1.shape[0]):
        jointProbs3_[i,:,:] = jointProbs1[i]*jointProbs2[:,:]
    arr = jointProbs3*(np.log2(jointProbs3)-np.log2(jointProbs3_))
    return np.mean(arr[np.isfinite(arr)]) 

# TCI for control
# mrsos_control[0].shape
# mrsos_std = np.zeros((mrsos_control[0].shape[1], mrsos_control[0].shape[2]))
# for i_lat in range(mrsos_control[0].shape[1]):
#     for j_lon in range(mrsos_control)
# https://stackoverflow.com/questions/58546999/calculate-correlation-in-xarray-with-missing-data
def linear_trend(x, y):
    #print(x.shape)
#     if np.sum(np.isnan(x))>0:
#         pf = np.empty((x.shape[0]))
#     else:
#         try:
    idx = np.isfinite(x) & np.isfinite(y)
    #print(x.shape)
    pf = np.polyfit(x[idx], y[idx], 1)
#         except:
#             pf = np.empty((x.shape[0]))
    return xr.DataArray(pf[0])

def compute_tci(sm, lh):
    print(sm.shape, lh.shape)
    sm_std = sm.std(dim='time')
    slopes = np.zeros_like(sm_std.values)

    x = sm
    y = lh
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)
    slopes     = cov/(xstd**2)
    intercept = ymean - xmean*slopes
    sm_std['slopes'] = (('lat', 'lon'), slopes)
    tci = sm_std.slopes*sm_std
    return tci
#mrsos_control[0].std(dim='time').plot()

def compute_aci(tas, hfss):
    sm = tas
    lh = hfss
    
    print(sm.shape, lh.shape)
    sm_std = sm.std(dim='time')
    slopes = np.zeros_like(sm_std.values)

    x = sm
    y = lh
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)
    slopes     = cov/(xstd**2)
    intercept = ymean - xmean*slopes
    sm_std['slopes'] = (('lat', 'lon'), slopes)
    tci = sm_std.slopes*sm_std
    return tci
#mrsos_control[0].std(dim='time').plot()

def soilm_memory(ds_soilm):
    threshold    = 1./np.exp(1.)
    soilm = ds_soilm.values
    smemory = np.zeros_like((ds_soilm.values[0,:,:]))
    ntim = 50 # ds_piClim_control_MPIESM_r1i1p1f1_mrso.mrso.values.shape[0]
    correlation_ = np.zeros((ntim, \
                        ds_soilm.shape[1], \
                        ds_soilm.shape[2]))    
    for tt in range(2,ntim):
        soilm_lagged  =  soilm[tt:,:,:]
        times = ds_soilm.time.values[tt:]
        lats = ds_soilm.lat.values
        lons = ds_soilm.lon.values
        ds = xr.Dataset({
        'soilm': xr.DataArray(
                    data   = soilm[:-tt],   # enter data here
                    dims   = ['time', 'lat', 'lon'],
                    coords = {'time': times, 'lat':lats, 'lon':lons},
                    ),
         'soilm_lagged': xr.DataArray(
                    data   = soilm_lagged,   # enter data here
                    dims   = ['time', 'lat', 'lon'],
                    coords = {'time': times, 'lat':lats, 'lon':lons},

                    )
                },
        )
        #print('lag = ', tt)
        x = ds['soilm']
        y = ds['soilm_lagged']
        n = y.notnull().sum(dim='time')
        xmean = x.mean(axis=0)
        ymean = y.mean(axis=0)
        xstd  = x.std(axis=0)
        ystd  = y.std(axis=0)

        #4. Compute covariance along time axis
        cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

        #5. Compute correlation along time axis
        cor   = cov/(xstd*ystd)
        correlation_[tt-2,:,:] = cor.values
    for i_lat in range(correlation_.shape[1]):
        for j_lon in range(correlation_.shape[2]):
            idx = np.where(correlation_[:,i_lat, j_lon] < 1/np.exp(1.))[0]
            #print(idx)
            #print(len(idx))
            #print(np.sum(np.isnan(soilm[:,i_lat,j_lon])))
            if len(idx)==0:
                smemory[i_lat, j_lon] = np.nan
            elif len(idx)==2:
                smemory[i_lat, j_lon] = np.nan
            else:
                print(idx)
                smemory[i_lat, j_lon] = idx[0]+1
    ds = xr.Dataset({
    'smemory': xr.DataArray(
                data   = smemory,   # enter data here
                dims   = [ 'lat', 'lon'],
                coords = {'lat':lats, 'lon':lons},
                )
            },
    )
            
    return ds
  
  def notaro_feedback_parameter(x,y):
    tau = 20 # daily data
    stptau = x.shift(time=tau)
    st = x
    atptau = y.shift(time=tau)
    n = atptau.notnull().sum(dim='time')
    xmean = st.mean(axis=0, skipna =True)
    ymean = atptau.mean(axis=0, skipna =True)
    xstd  = st.std(axis=0, skipna =True)
    ystd  = atptau.std(axis=0, skipna =True)
    cov_num = np.sum((st - xmean)*(atptau - ymean), axis=0)/(n) 

    n = stptau.notnull().sum(dim='time')
    xmean = st.mean(axis=0, skipna =True)
    ymean = stptau.mean(axis=0, skipna =True)
    xstd  = st.std(axis=0, skipna =True)
    ystd  = stptau.std(axis=0, skipna =True)
    cov_den = np.sum((st - xmean)*(stptau - ymean), axis=0)/(n) 

    nfp = cov_num/cov_den
    
    return nfp
  
  def zengs_gamma(x,y):
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)

    #4. Compute covariance along time axis
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

    #5. Compute correlation along time axis
    cor   = cov/(xstd*ystd)
    
    cor = cor * (xstd/ystd)
    return cor
  
  
