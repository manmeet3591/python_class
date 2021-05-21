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
