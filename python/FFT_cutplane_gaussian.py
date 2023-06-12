import numpy as np
import sys
import os
import matplotlib.pyplot as plt 
from scipy.fftpack import fft, ifft
import matplotlib.animation as animation

rho0=1.2
T0=300
M=0.4
U0=M*np.sqrt(1.4*287*300)
P0= 1/1.4/M/M

#sys.path.append('/home/student.unimelb.edu.au/shubham1/Dropbox/PhDWork/python_codes/slicedatareadWriteHipstar')
path = "/data/cephfs/punim0895/hwu1/LES_new_SS/M04/Subspaces/"
import lib_read_HiPSTAR as hp
import os
blockid =2
blockid2 =3
blockname = 'xy_cutplane'
blockname2 = 'temp_surf'
nvar=5
i_index=0
j_index=0
k_index=0
istart = i_index
iend   = i_index
jstart = j_index
jend   = j_index
kstart = k_index
kend   = k_index
SUBVOL=[  [[blockid],[istart,iend],[jstart,jend],[None]] ]
SUBVOL=[]
grid = hp.get_grids(wdir=path,subvol=SUBVOL,read_subspaces=True,subspaces_list=[[blockid],[blockname]])
grid2= hp.get_grids(wdir=path,subvol=SUBVOL,read_subspaces=True,subspaces_list=[[blockid2],[blockname]])
iteration=range(531000,800000,1000)  
iteration=range(1915200,3175000,200)  
iteration=range(1800100,1800100+149900,100)   # M02
# iteration=range(970100,970100+299900,100)  # M03
iteration=range(790100,790100+49900,100) # M04
#iteration=range(710100,710100+149900,100)   # M05
DATA=np.zeros((grid[0].idim,grid[0].jdim,grid[0].kdim,nvar,len(iteration)))
DATA2=np.zeros((grid2[0].idim,grid2[0].jdim,grid2[0].kdim,nvar,len(iteration)))
################      READING AND WRITING SUBVOLUMES        ##############################################
for time in range(len(iteration)):
#    print(time)
    flow = hp.get_data(grid,wdir=path,tstep=iteration[time],ss_var=nvar)
    flow2 = hp.get_data(grid2,wdir=path,tstep=iteration[time],ss_var=nvar)
    DATA[:,:,:,:,time]=flow[0].data
    DATA2[:,:,:,:,time]=flow2[0].data
######################processing##################################################################################
Nt = len(iteration)
rho = DATA[:,:,0,0,:]
T   = DATA[:,:,0,4,:]
# rho = DATA[:,:,0,0,:Nt] + DATA[:,:,0,0,Nt:2*Nt]
# T =   DATA[:,:,0,4,:Nt] + DATA[:,:,0,4,Nt:2*Nt]
p   = rho*T*P0 
rho2 = DATA2[:,:,0,0,:]
T2   = DATA2[:,:,0,4,:]
# rho2 = DATA2[:,:,0,0,:Nt] + DATA2[:,:,0,0,Nt:2*Nt]
# T2 =   DATA2[:,:,0,4,:Nt] + DATA2[:,:,0,4,Nt:2*Nt]
p2   = rho2*T2*P0 
pmean_temp = np.mean(p,axis=-1)
pmean = pmean_temp[...,np.newaxis]
pfluc = np.subtract(p,pmean)
pmean2_temp = np.mean(p2,axis=-1)
pmean2 = pmean2_temp[...,np.newaxis]
pfluc2 = np.subtract(p2,pmean2)
#########################################################################################################
npx = grid[0].npx
npy = grid[0].npy
npx2 = grid2[0].npx
npy2 = grid2[0].npy
print npx, npy, npx2, npy2
x = grid[0].data[:,:npy-2,0,0]
y = grid[0].data[:,:npy-2,0,1]
x2 = grid2[0].data[:,:npy2-2,0,0]
y2 = grid2[0].data[:,:npy2-2,0,1]
pfluc_mod = pfluc[:,:npy-2,:]
pfluc2_mod = pfluc2[:,:npy2-2,:]
print np.argwhere(np.isnan(pfluc2_mod[:,:,1]))
print np.argwhere(np.isnan(pfluc_mod[:,:,1]))

levels=np.linspace(np.minimum(np.min(pfluc_mod[:,:,4]),np.min(pfluc2_mod[:,:,4])),
                   np.maximum(np.max(pfluc_mod[:,:,4]),np.max(pfluc2_mod[:,:,4])),7)
print levels
levels=np.linspace(-0.3,0.1,10)
levels=None
fig = plt.figure(figsize=(20,14))    
plt.contourf(x,y,pfluc_mod[:,:,4],256,levels=levels,cmap='coolwarm')
# plt.contourf(x2,y2,pfluc2_mod[:,:,4],256,levels=levels,cmap='coolwarm')
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.show()
########################################################################################################
gridSize = grid[0].npx*(grid[0].npy-2)
pfft = np.zeros((npx,npy-2,Nt),dtype=complex)
pfft2 = np.zeros((npx2,npy2-2,Nt),dtype=complex)

pfft =  fft(pfluc_mod , axis = -1)
pfft2 = fft(pfluc2_mod, axis = -1)
#####################################################################################################
dtnon = 5e-5
Lref = 100000*1.5e-5/U0
tconv = Lref/U0
dt = dtnon*tconv*100
print dt
tf = np.linspace(0, 1.0/(1.0*dt), Nt) *Lref/U0
print Nt
fig = plt.figure(figsize=(6.5,5.5))    
WPS=(2*np.abs(pfft[500,0,1:Nt//2]))
cba1 = plt.loglog(tf[1:Nt//2],WPS)
plt.tight_layout()
plt.show()
########################################################################################################
#####  freq selction matrix
strouhal_matrix=[0.5,2 2,8,8,11,0.5,11,11,20,11,30,20,50,30,100,11,100]
len_strouhal =len(strouhal_matrix)/2
for freq in range(len_strouhal):
    freq1 = strouhal_matrix[freq*2]   # selction fo frequencies in pair
    freq2 = strouhal_matrix[freq*2+1]
    dirpath = 'M'+ str(M) + '_' + str(freq1) + '-' + str(freq2)
    os.mkdir(dirpath)
# selction of frequency range(7-10Hz) for time domain conversion
    index=np.where(tf>freq1)         #   index with frequency lower than this will be removed
    index1=np.where(tf<freq2)         #   index with frequency lower than this will be removed
    print index[0][0]
    print index1[0][-1]
############################ GAUSSIAN WINDOW DEFINITION ##############################################
# Introduce Gaussian function
    sigma=5
    window=np.zeros(len(tf))
    window[0:index[0][0]] = np.exp(-(tf[0:index[0][0]]-freq1)**2/(2*sigma**2))
    window[index[0][0]:index1[0][-1]] = 1
    window[index1[0][-1]:Nt//2]  = np.exp(-(tf[index1[0][-1]:Nt//2]-freq2)**2/(2*sigma**2))
    window[Nt//2:-index1[0][-1]] = np.exp(-(tf[Nt//2:-index1[0][-1]]+freq2-tf[-1])**2/(2*sigma**2))
    window[-index1[0][-1]:-index[0][0]] = 1
    window[-index[0][0]:] = np.exp(-(tf[-index[0][0]:]+freq1-tf[-1])**2/(2*sigma**2))
    plt.figure()
    plt.plot(tf,window)
    plt.title(r"Gaussian window ($\sigma$=7)")
    plt.ylabel("Amplitude")
    plt.xlabel("Sample")

############################### INVERSE FFT #############################################
# IFFT for BLOCK 2
    pifft = np.zeros((npx,npy-2,Nt),dtype=complex)   # converted to time-domain
    pifftabs = np.zeros((npx,npy-2,Nt),dtype=complex)   # extract the abs values for one-sided spectrum
    pfft_gauss=pfft*window;            # multiplied to gaussian function

    pifft = ifft(pfft_gauss, axis = -1)     #  Bring the data back into time-domain
    pifftreal=np.real(pifft);  # extract only the real data
#######################################################################################################
# IFFT for BLOCK 3
    pifft2 = np.zeros((npx2,npy2-2,Nt),dtype=complex)
    pifftabs2 = np.zeros((npx2,npy2-2,Nt),dtype=complex)
    pfft2_gauss=pfft2*window;            # multiplied to gaussian function

    pifft2 = ifft(pfft2_gauss, axis = -1)     #  Bring the data back into time-domain
    pifftreal2=np.real(pifft2);  # extract only the real data
#######################################################################################################

# For specific frequency
# levels=np.linspace((np.min(pifftreal)),(np.max(pifftreal)),256).real
    levels=np.linspace(-0.14,0.07,256)
    levels=np.linspace(-1,1,256)
    factor=5e-3

# levels=None
    for j in range(200):
        i=j*2
        fig,ax = plt.subplots(figsize=(30,15))
        pcm=plt.contourf(x,y,pifftreal[:,:,i]/factor,256,levels=levels,cmap='gray',extend='both')
        pcm=plt.contourf(x2[:,:300],y2[:,:300],pifftreal2[:,:300,i]/factor,256,levels=levels,cmap='gray',extend='both')
        if i<10:
            dtemp1 = '000' + str(i)
        elif i<100:
            dtemp1 = '00' + str(i)
        elif i<1000:
            dtemp1 = '0' + str(i)
        else:
            dtemp1 = str(i)
            
        dtprint = i*dt
        dtemp=np.format_float_positional(np.float16(dtprint), unique=False, precision=7)
        plt.title('t='+ dtemp)
        cbar=plt.colorbar(pcm)
        plt.savefig('./'+ dirpath +'/M' +str(M)+'_.'+dtemp1+'.png',dpi=300,bbox_inches="tight")
        plt.close(fig)
