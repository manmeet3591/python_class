arr1 = ds_tp_era5_[0].copy()
arr2 = ds_dl_[0].copy()
thres = 0.5
def csi_ets(arr1, arr2, thres):
    
    arr1_  = arr1>thres
    arr2_  = arr2>thres
    joint_bool = (arr1>thres)*(arr2>thres)
    hits = np.sum(joint_bool)
    arr3_ = arr1_
    arr3_[(arr1_==True) & (arr2_==True)] = False
    misses = np.sum(arr3_*arr1_)
    arr4_ = arr2_
    arr4_[(arr1_==True) & (arr2_==True)] = False
    false_alarms = np.sum(arr4_*arr2_)

    csi = hits / (hits+misses+false_alarms)
    e = arr1.shape[0]
    ets = (hits - e)/(hits+misses+false_alarms-e)

    return csi, ets
csi, ets = csi_ets(arr1, arr2, thres)
