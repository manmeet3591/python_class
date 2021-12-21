from collections import defaultdict
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
        res+=v*np.log2(v)
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
    hx, bx= np.histogram(x, bins=int(np.sqrt(x.size/5)),density=True)

    Py= compute_distribution(y)
    Px= compute_distribution(np.digitize(x,bx))

    res= 0
    for ey in set(y):
        # P(X | Y)
        x1= x[y==ey]
        condPxy= compute_distribution(np.digitize(x1,bx))

        for k, v in condPxy.items():
            res+= (v*Py[ey]*(np.log2(Px[k]) - np.log2(v*Py[ey])))
    return res
            
def mutual_information(x,y):
    return entropy(y) - conditional_entropy(x,y)
    
def mutual_information_ufunc(x, y, dim='time'):
    return xr.apply_ufunc(
        mutual_information,
        x,
        y,
        input_core_dims=[[dim], [dim]],
        dask="parallelized",
        output_dtypes=[float],
    )
  
 mi = mutual_information_ufunc(x.chunk({'time': -1}),y.chunk({'time': -1}))
