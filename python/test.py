
# coding: utf-8

# In[1]:


# General modules and packages
import matplotlib.pyplot as plt
import numpy as np
import shutil
import os
import pylab as pl
#pl.rc('text', usetex=True)
#%config InlineBackend.figure_format = 'svg'
#%matplotlib inline 

# other downloaded modules
#from matplotlib2tikz import save as tikz_save

#own modules
#import tools
#import lib_hipstar as lh
#import lines
#import deriv
from scipy.interpolate import interp1d
import math as mt
import multiprocessing as mp


# ## Calculate the shift required in the data

# In[2]:


dt_simulation =   0.00004188790204765447 #0.0001745329251985603 #0.00004188790204765447
save_interval =   100
Omega         =   0.666666667
modes         =   256
t_phys        =   dt_simulation*save_interval
theta_travel  =   Omega*t_phys
theta_travel_rad = (180./(np.pi))*theta_travel
print theta_travel_rad
grid_point_travel= ((2.*modes+2)/360.)*theta_travel_rad
print grid_point_travel
gpt=int(grid_point_travel)
print gpt


# In[3]:


print theta_travel_rad
Theta_Rad=((np.pi/180.)*theta_travel_rad)
print Theta_Rad


# In[4]:


#wdir = '/home/sainideepak/Phd_work/JFM_paper/Magnus/Ra1e09/Ro_0.75_test/FLOW/'
wdir = '/home/sainideepak/Phd_work/JFM_paper/Magnus/Ra1e09/Ro_0.75_test/FLOW/FLOW_d_new/'

# In[5]:


filename = wdir+'FLOW_phys_1_'


# In[6]:

ts_start = 5040011
ts_end =   5044911
stepsize = 100
files_no = (ts_end - ts_start)/stepsize + 1
dt=0.000314159265359 #0.000314159265359 #0.000010030627885024883
samplingrate = 1.0/(stepsize*dt)
data_type2 = np.dtype('float32')
nxp = 380
nyp = 380
nzp = 2*256+2
x     = np.zeros((nxp,nyp,nzp))
y     = np.zeros((nxp,nyp,nzp))
z     = np.zeros((nxp,nyp,nzp))
theta     = np.zeros((nxp,nyp,nzp))
print files_no



# In[7]:

print('Reading Grid')

f = open(wdir+'FLOW_phys_GRID_1.xyz', 'rb')
data_type1 = np.dtype('int32')
nxp,nyp,nzp = np.fromfile(f,dtype=data_type1,count=3)
print nxp,nyp,nzp
for k in range(nzp):
    for j in range(nyp):
        x[:,j,k] = np.fromfile(f,dtype=data_type2,count=nxp)
for k in range(nzp):
    for j in range(nyp):
        y[:,j,k] = np.fromfile(f,dtype=data_type2,count=nxp)        
for k in range(nzp):
    for j in range(nyp):
        z[:,j,k] = np.fromfile(f,dtype=data_type2,count=nxp)
f.close()


# In[8]:
cpu=8
pool=mp.Pool()

def read_data2(files_no):
    print(files_no)


def read_data(files_no):
    data_type2 = np.dtype('float32')
    T_midx =  np.empty((nyp,nzp),dtype=data_type2)
    T =  np.empty((nxp,nyp,nzp),dtype=data_type2)
    istep = ts_start+(files_no)*stepsize
    f = open(filename+str(istep)+'.raw','rb')
    f.seek(28+4*nxp*nyp*nzp*4,0)
    for k in range(nzp):
        for j in range(nyp):
            T[:,j,k] = np.fromfile(f,dtype=data_type2,count=nxp)
    f.close()
    T_midx[:,:]=T[nxp/2,:,:]
    return T_midx
                                                                            


print files_no
results=pool.map(read_data2,range(files_no))
pool.close()
pool.join()

