def pattern_cor(ds1, ds2):
    # Assume latxlon grid
    nlat = ds1.lat.values.shape[0]
    nlon = ds1.lon.values.shape[0]
    weight = np.zeros((nlat, nlon))
    for i_lat in range(nlat):
        weight[:,0] = np.cos(ds1.lat.values[i_lat]*(np.pi/180))
    ds1_ = np.multiply(ds1.values[0,:,:],weight)
    ds1_sum = np.sum(ds1_.flatten())
    ds1_mean = ds1_sum/(np.sum(weight.flatten()))

    ds2_ = np.multiply(ds2.values[0,:,:],weight)
    ds2_sum = np.sum(ds2_.flatten())
    ds2_mean = ds2_sum/(np.sum(weight.flatten()))
    
    xycov = 0.0
    xanom2 = 0.0
    yanom2 = 0.0
    w = weight
    x = ds1.values[0,:,:]
    y = ds2.values[0,:,:]
    xave = ds1_mean                 
    yave = ds2_mean
    print(weight.shape, x.shape, y.shape, xave, yave)             
    for ml in range(nlat):
        for nl in range(nlon):
            xycov  = xycov  + w[ml,nl]*(x[ml,nl]-xave)*(y[ml,nl]-yave)
            xanom2 = xanom2 + w[ml,nl]*(x[ml,nl]-xave)**2
            yanom2 = yanom2 + w[ml,nl]*(y[ml,nl]-yave)**2
                      
    r   = xycov/(np.sqrt(xanom2)*np.sqrt(yanom2))
    return r
