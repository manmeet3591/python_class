#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras.layers import *
import tensorflow.keras.backend as K
import seaborn as sns
import pickle
#from src.score import * # load_test_data defined here
from collections import OrderedDict

from tensorflow.keras.layers import Input, UpSampling3D, AveragePooling3D, concatenate, ReLU, Reshape, Concatenate,     Permute
#from DLWP.custom import CubeSpherePadding2D, CubeSphereConv2D
from tensorflow.keras.layers import Reshape, Concatenate, Permute
from tensorflow.keras.models import Model
import glob
from tensorflow.keras.models import load_model

import os
import numpy as np

from tensorflow.keras.models import Sequential, model_from_json
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Activation
import cv2
from patchify import * 
#import tensorly as tl


# In[2]:


def mse(y, t):
    return np.mean(np.square(y - t))

def psnr(y, t, max_=200.0):
    return 20 * np.log10(max_) - 10 * np.log10(mse(y, t))

def ssim(x, y, max_=200.0):
    mu_x = np.mean(x)
    mu_y = np.mean(y)
    var_x = np.var(x)
    var_y = np.var(y)
    cov = np.mean((x - mu_x) * (y - mu_y))
    c1 = np.square(0.01 * max_)
    c2 = np.square(0.03 * max_)
    return ((2 * mu_x * mu_y + c1) * (2 * cov + c2)) / ((mu_x**2 + mu_y**2 + c1) * (var_x + var_y + c2))


# In[3]:


data_dir='/lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/austin/'
file_ = glob.glob(data_dir+'GSMaP.v7_0.10deg_????-*_austin.nc')


# In[4]:


ds = xr.open_mfdataset(file_, combine='by_coords')


# In[5]:


ds_test  = ds.sel(time=slice('2010', '2020'))
min_ = 0.0 # min_max.values[0]
max_ = 200.0 #min_max.values[1]


# In[6]:


weight_filename = 'srcnn_weight_2008.hdf5'
model = Sequential()
model.add(Conv2D(64,9,padding='same',input_shape=(20,20,1)))
model.add(Activation('relu'))
model.add(Conv2D(32,1,padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(1,5,padding='same'))
# optimizer = Adam(lr=0.001)
# model.compile(optimizer=optimizer, loss='mean_squared_error', metrics=['accuracy'])
model.load_weights(os.path.join('./model/',weight_filename))


# In[7]:


ds_test  = ds.sel(time=slice('2010', '2020'))
year = 2010
data_ = ds_test.pr.sel(time=slice(str(year), str(year))).values


# In[8]:


data_n = 2*((data_[:,:,:]-min_)/(max_-min_))
print(data_.shape, data_n.shape)

X = []
y = []
for i_ in range(data_n.shape[0]):
    img = data_n[i_,:,:]
    if True:#np.sum(img)>0:
        scale = 2
        res = cv2.resize(img, dsize=(img.shape[0]*scale, img.shape[1]*scale), interpolation=cv2.INTER_CUBIC)
        res_ = cv2.resize(res, dsize=(int(res.shape[0]*(1/scale)),                                           int(res.shape[1]*(1/scale))), interpolation=cv2.INTER_CUBIC)
            #data_n_[i_,:,:] = res_
        X.append(img)
        y.append(res_)
X = np.asarray(X)
y = np.asarray(y) 


# In[134]:


sub_images_ = patchify(X[0], (20,20), step=10)
sub_labels_ = patchify(X[0], (20,20), step=10)
fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(11.69,8.27))
ax[0].imshow(sub_images_[0,0])
sub_pred1 = model.predict(sub_images_[0,0,np.newaxis,:,:,np.newaxis])[0]
sub_pred1[sub_pred1<=1e-4] = 0.0
ax[1].imshow(sub_pred1)
print(np.sum(sub_images_[0,0]))
print(np.sum(model.predict(sub_images_[0,0,np.newaxis,:,:,np.newaxis])[0]))


# In[110]:


fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(11.69,8.27))
ax[0].imshow(sub_images_[0,1])
sub_pred2 = model.predict(sub_images_[0,1,np.newaxis,:,:,np.newaxis])[0]
sub_pred2[sub_pred2<=1e-4] = 0.0
ax[1].imshow(sub_pred2)
print(np.sum(sub_images_[0,1]))
print(np.sum(sub_labels_[0,1]))
print(np.sum(model.predict(sub_images_[0,1,np.newaxis,:,:,np.newaxis])[0]))
print(np.min(sub_pred2), np.max(sub_pred2), np.min(sub_images_[0,1]), np.max(sub_images_[0,1]))


# In[105]:


fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(11.69,8.27))
ax[0].imshow(sub_images_[1,0])
sub_pred3 = model.predict(sub_images_[1,0,np.newaxis,:,:,np.newaxis])[0]
sub_pred3[sub_pred3<=1e-4] = 0.0
ax[1].imshow(sub_pred3)
print(np.sum(sub_images_[1,0]))
print(np.sum(sub_labels_[1,0]))
print(np.sum(model.predict(sub_images_[1,0,np.newaxis,:,:,np.newaxis])[0]))
print(np.min(sub_pred3), np.max(sub_pred3), np.min(sub_images_[1,0]), np.max(sub_images_[1,0]))


# In[139]:


fig,ax = plt.subplots(ncols=2,nrows=1, figsize=(11.69,8.27))
ax[0].imshow(sub_images_[1,1])
sub_pred4 = model.predict(sub_images_[1,1,np.newaxis,:,:,np.newaxis])[0]
sub_pred4[sub_pred4<=1e-4] = 0.0
ax[1].imshow(sub_pred4)
print(np.sum(sub_images_[1,1]))
print(np.sum(sub_labels_[1,1]))
print(np.sum(model.predict(sub_images_[1,1,np.newaxis,:,:,np.newaxis])[0]))
print(np.min(sub_pred4), np.max(sub_pred4), np.min(sub_images_[1,1]), np.max(sub_images_[1,1]))
print(sub_pred4.shape)


# In[159]:


sub_images_ = patchify(X[0], (20,20), step=10)
patches =np.copy(sub_images_)
sub_pred = model.predict(sub_images_[0,0,np.newaxis,:,:,np.newaxis])[0]
patches[0,0] = sub_pred[:,:,0]

sub_pred = model.predict(sub_images_[0,1,np.newaxis,:,:,np.newaxis])[0]
patches[0,1] = sub_pred[:,:,0]

sub_pred = model.predict(sub_images_[1,0,np.newaxis,:,:,np.newaxis])[0]
patches[1,0] = sub_pred[:,:,0]

sub_pred = model.predict(sub_images_[1,1,np.newaxis,:,:,np.newaxis])[0]
patches[1,1] = sub_pred[:,:,0]
#patches[sub_images_<=1e-5]=0.0
patches[sub_images_==0.0] = 0.0
reconstructed_image = unpatchify(patches, image.shape)

#print(np.min(patches[0,0]),      


# In[160]:


plt.imshow(reconstructed_image)


# In[161]:


print(psnr(reconstructed_image, y[0]), psnr(X[0], y[0]))


# In[144]:


plt.imshow(sub_images_[0,1])


# In[145]:


#print(X[0].shape)
X_pred = np.copy(X)
for i_x in range(X.shape[0]):
    for i in range(0,X[i_x].shape[0]-10,10):
        for j in range(0,X[i_x].shape[1]-10,10):
            print(i, i+20,j, j+20)
            if np.sum(X[i_x:i_x+1,i:i+20,,,,,,,j:j+20])==0.0:
                continue
            X_pred[i_x,i:i+20, j:j+20] = model.predict(X[i_x:i_x+1,i:i+20,j:j+20,np.newaxis])[:,:,0]


# In[ ]:


psnr(y, X_pred), psnr(X,y)


# In[ ]:





# In[96]:


# sub_pred4[sub_pred4<=1e-4] = 0.0
# if np.sum(sub_images_)
X_pred = np.copy(X)
sub_pred_list = []
# Transform to sub images and predict
for i_ in range(1):#X.shape[0]):
    print(i_)
    image = X[i_]
    sub_images_ = patchify(image, (20,20), step=10)

    image = y[i_]
    sub_labels_ = patchify(image, (20,20), step=10)
    sub_pred_ = np.zeros_like((sub_images_))
    print(sub_images_.shape)
    #print(sub_images_.shape)
    for j_ in range(sub_images_.shape[0]):
        for k_ in range(sub_images_.shape[1]):
            print(j_,k_)
            #if np.sum(sub_images_[j_,k_])>0.0:
            sub_pred_[j_,k_,:,:] = model.predict(sub_images_[j_,k_,np.newaxis,:,:,np.newaxis])[:,:,0]
            sub_pred_list.append(sub_pred_[j_,k_,:,:])
        #sub_pred_[sub_pred_<=1e-4] = 0.0
            #sub_pred_[j_,k_,:,:][sub_pred_[j_,k_,:,:]<=1e-5]=0.0
            #sub_pred_[j_,k_,:,:][sub_images_[j_,k_,:,:]==0.0]=0.0
        
            #print(sub_pred_.shape)
    pred_reconstructed = unpatchify(sub_pred_, image.shape)
    print(pred_reconstructed.shape)
    X_pred[i_] = pred_reconstructed


# In[99]:


plt.imshow(sub_pred_list[4])


# In[98]:


import matplotlib.pyplot as plt
fig,ax = plt.subplots(ncols=3,nrows=1, figsize=(11.69,8.27))
data_ = [y[0], X[0], X_pred[0]]
tits_ = ['Ground Truth', 'Cubic', 'SRCNN']
for i in range(3):
    ax[i].imshow(data_[i])
    ax[i].set_title(tits_[i])


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




