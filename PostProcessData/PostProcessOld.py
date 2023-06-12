# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:09:01 2016
This file contains functions to post-process data from an OpenFoam field, which are generated from
the sample command. 

@author: harshal
"""
#Global variable
gamma = 1.4
import copy
#class line(object):
import numpy as np


def mixedOutQuantity(rho,p,u,v,y):
    "This function returns mixed out quantities. "
    import numpy as np
    from math import sqrt,atan
   
    gam1=gamma-1.0
    mdot=0.0
    m1_mix=0.0
    m2_mix=0.0
    energ=0.0
    l = np.abs(y[0]-y[-1])
    
    for j in range(1,len(y)):
        Delta = y[j-1]-y[j]
        mdot=mdot+((rho[j]*u[j]+rho[j-1]*u[j-1])*Delta/2)
        m1_mix=m1_mix+((rho[j]*u[j]*u[j]+p[j]+ rho[j-1]*u[j-1]*u[j-1]+p[j-1])*Delta/2)
        m2_mix=m2_mix+((rho[j]*u[j]*v[j]+ rho[j-1]*u[j-1]*v[j-1])*Delta/2)
        energ=energ+((u[j]*(gamma/gam1*p[j]+(u[j]*u[j]+ v[j]*v[j])*rho[j]*0.5) + u[j-1]*(gamma/gam1*p[j-1]+(u[j-1]*u[j-1]+ v[j-1]*v[j-1])*rho[j-1]*0.5))*Delta/2)
      
    m0_store=mdot/l
    m1_mix_out=m1_mix/l
    m2_mix_out=m2_mix/l
    ener_mix_out=energ/l
    
    neg_p_half = m1_mix_out/m0_store*gamma/(gamma+1)
    neg_q = 2./m0_store*gam1/(gamma+1.)*(m2_mix_out**2/m0_store*0.5-ener_mix_out)
    u_mix_out = neg_p_half-sqrt(neg_p_half**2+neg_q)
    rho_mix_out = m0_store/u_mix_out
    v_mix_out=m2_mix_out/rho_mix_out/u_mix_out
    p_mix_out=m1_mix_out-rho_mix_out*u_mix_out*u_mix_out
    
    t_mix_out = p_mix_out/(rho_mix_out*287)
    c_mix_out = (u_mix_out**2+v_mix_out**2)**0.5
    ma_mix_out = c_mix_out / (gamma*p_mix_out/rho_mix_out)**0.5
    pt_mix_out = p_mix_out*(1+0.2*ma_mix_out**2)**(1.4/0.4)
    tt_mix_out = t_mix_out*(1+0.2*ma_mix_out**2)
    mu_mix_out = t_mix_out**1.5*1.458*(10**-6)/(t_mix_out+110.4)
    angle_mix_out = atan(v_mix_out/u_mix_out)*45/atan(1.)
    
    return [p_mix_out, pt_mix_out, rho_mix_out, u_mix_out, v_mix_out, c_mix_out, ma_mix_out, t_mix_out, tt_mix_out, mu_mix_out, angle_mix_out]
    
def xtractData(directory,time1,location,model,purpose):
    "This function is used to extract the mixed out losses."
    "It also extracts various field quantities in the vertical direction (in, out, up, down)"
    import numpy as np
    fileVelo ='a'
    fileOther = 'b'
    if model == 'SA':
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_nuTilda_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_nuTilda_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_nuTilda_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_U.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_nuTilda_rho_T.xy'
        else:
            print 'Please enter correct location'
                 
    elif model == 'SST':
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_omega_k_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_omega_k_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_omega_k_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_U.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_omega_k_rho_T.xy'
        else:
            print 'Please enter correct location'
            
    elif model == 'LKE':
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_omega_kt_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_omega_kt_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_omega_kt_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_U.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_omega_kt_rho_T.xy'
        else:
            print 'Please enter correct location'
        
    else:  #The model is kEpsilon
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_epsilon_k_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_epsilon_k_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_epsilon_k_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_U.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_epsilon_k_rho_T.xy'
        else:
            print 'Please enter correct location'  
        
    velocity = np.loadtxt(fileVelo)
    other    = np.loadtxt(fileOther) # better to put these above.  and return proper quanties - instead of other - causes confusion. 
     
#    0 - y; 1 p; 2 -nuTilda; 3 rho; 4 T; - SA Model 
#   0 - y; 1 p; 2 - omega / epsilon; 3 k; 4 rho; 5 T; - other 
    MagVelo = (velocity[:,1]**2 + velocity[:,2]**2)**0.5 
    
    
    if purpose == 'MixedOut':
        if model == 'SA':
            mixedOutQuant = mixedOutQuantity(other[:,3],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        elif model =='LKE':
            mixedOutQuant = mixedOutQuantity(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        else:
            mixedOutQuant = mixedOutQuantity(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        return mixedOutQuant
    elif purpose == 'wakeLoss':
        if model == 'SA':
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,3])**0.5
            PStag   = other[:,1]*(1+0.2*MachNo**2)**(1.4/0.4)
        elif model == 'LKE':
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,4])**0.5
            PStag   = other[:,1]*(1+0.2*MachNo**2)**(1.4/0.4) 
        else:
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,4])**0.5
            PStag   = other[:,1]*(1+0.2*MachNo**2)**(1.4/0.4)
        return [other[:,0], other[:,1],PStag,MagVelo,velocity[:,1],velocity[:,2]]
    
    else:  #Time convergence
        return [other[:,0], other[:,1], MagVelo,other[:,2],other[:,3],other[:,4],other[:,5]]
        
        
        
def xtractDataTimeAv(directory,time1,location,model,purpose):
    "This function is used to extract the mixed out losses for time averaged flow fields."
    "It also extracts various field quantities in the vertical direction (in, out, up, down)"
    import numpy as np
    fileVelo ='a'
    fileOther = 'b'
    if model == 'SA':
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_pMean_nuTildaMean_rhoMean_TMean.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_pMean_nuTildaMean_rhoMean_TMean.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_pMean_nuTildaMean_rhoMean_TMean.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_UMean.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_pMean_nuTildaMean_rhoMean_TMean.xy'
        else:
            print 'Please enter correct location'
                 
    elif model == 'SST':
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_omega_k_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_omega_k_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_omega_k_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_UMean.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_omega_k_rho_T.xy'
        else:
            print 'Please enter correct location'


    elif model == 'LKE':
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_pMean_omegaMean_kt_rho_TMean.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_pMean_omegaMean_kt_rho_TMean.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_pMean_omegaMean_kt_rho_TMean.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_UMean.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_pMean_omegaMean_kt_rho_TMean.xy'
        else:
            print 'Please enter correct location'        
    else:  #The model is kEpsilon
        if location =='in':   
            fileVelo  = directory+'sets/'+str(time1)+'/inlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_epsilon_k_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_epsilon_k_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_UMean.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_epsilon_k_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_UMean.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_epsilon_k_rho_T.xy'
        else:
            print 'Please enter correct location'  
        
    velocity = np.loadtxt(fileVelo)
    other    = np.loadtxt(fileOther) # better to put these above.  and return proper quanties - instead of other - causes confusion. 
     
#    0 - y; 1 p; 2 -nuTilda; 3 rho; 4 T; - SA Model 
#   0 - y; 1 p; 2 - omega / epsilon; 3 k; 4 rho; 5 T; - other 
    MagVelo = (velocity[:,1]**2 + velocity[:,2]**2)**0.5 
    
    
    if purpose == 'MixedOut':
        if model == 'SA':
            mixedOutQuant = mixedOutQuantity(other[:,3],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        else:
            mixedOutQuant = mixedOutQuantity(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        return mixedOutQuant
    elif purpose == 'wakeLoss':
        if model == 'SA':
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,3])**0.5
            PStag   = other[:,1]*(1+0.2*MachNo**2)**(1.4/0.4)
            
        else:
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,4])**0.5
            PStag   = other[:,1]*(1+0.2*MachNo**2)**(1.4/0.4)
        return [other[:,0], other[:,1],PStag,MagVelo,velocity[:,1],velocity[:,2]]
    
    else:  #Time convergence
        return [other[:,0], other[:,1], MagVelo,other[:,2],other[:,3],other[:,4],other[:,5]]
    
#    % Plots the Cp for given pressure profiles over the blade. 
def plotCp(directory,time1,mark,n,label2,model,fine='no'):  
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(n)  
    filename = directory+'surfaces/'+str(time1)+'/p_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0]
    p = A[:,3]
    print len(p)
    mixedOutUp = xtractData(directory,time1,'up',model,'MixedOut')
    mixedOutDown = xtractData(directory,time1,'down',model,'MixedOut')  
    if fine=='yes':
        p2 = np.append(p[0:287], np.append(p[646:1291], np.append(p[288:645],p[0])))
        p1 = np.append(p2[81:877],np.append(p2[878:1291],p2[0:82]))
        x2 = np.append(x[0:287], np.append(x[646:1291], np.append(x[288:645],x[0])))
        x1 = np.append(x2[81:877],np.append(x2[878:1291],x2[0:82]))
    else:
        p2 = np.append(p[0:190], np.append(p[430:859], np.append(p[191:429],p[0]))) #reorder the data around the blade. 
        p1 = np.append(p2[50:585],np.append(p2[586:858],p2[0:49]))
        x2 = np.append(x[0:190], np.append(x[430:859], np.append(x[191:429],x[0])))
        x1 = np.append(x2[50:585],np.append(x2[586:858],x2[0:49])) #data is now reorder so that 0-585 is suction side, 586-858  is pressure side
    cp = (p1- mixedOutDown[0]) / (mixedOutUp[1]-mixedOutDown[0])
    plt.plot(x1/0.0859693,cp,mark,linewidth=2,label = label2)
    plt.legend()
    plt.xlim(-0.02,1)
    plt.ylim(-0.7,1.01)
#    plt.title('Loading on the T106A blade for different inlet angles')
    plt.xlabel('x/C$_{ax}$',fontsize=20)
    plt.ylabel('C$_p$',fontsize=20)
    
def plotyPlus(directory,time1,mark,n,label2,fine='no'):  
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(n)  
    filename = directory+'surfaces/'+str(time1)+'/yPlus_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0]
    yPlus = A[:,3]
    print yPlus
    ss = mark
    ps = mark + '--'
    if fine=='yes':
        yPlus2 = np.append(yPlus[0:287], np.append(yPlus[646:1291], np.append(yPlus[288:645],yPlus[0]))) #reorder the data around the blade. 
        yPlus1 = np.append(yPlus2[81:877],np.append(yPlus2[878:1291],yPlus2[0:80]))
        x2 = np.append(x[0:287], np.append(x[646:1291], np.append(x[288:645],x[0])))
        x1 = np.append(x2[81:877],np.append(x2[878:1291],x2[0:80]))
        plt.plot(x1[:810]/0.0859693,yPlus1[:810],ss,linewidth=2,label = label2 + ' - Suction Side')
        plt.plot(x1[811:]/0.0859693,yPlus1[811:],ps,linewidth=2,label = label2 + ' - Pressure Side')
    else:            
        yPlus2 = np.append(yPlus[0:190], np.append(yPlus[430:859], np.append(yPlus[191:429],yPlus[0]))) #reorder the data around the blade. 
        yPlus1 = np.append(yPlus2[50:585],np.append(yPlus2[586:858],yPlus2[0:49]))
        x2 = np.append(x[0:190], np.append(x[430:859], np.append(x[191:429],x[0])))
        x1 = np.append(x2[50:585],np.append(x2[586:858],x2[0:49]))
        plt.plot(x1[:540]/0.0859693,yPlus1[:540],ss,linewidth=2,label = label2 + ' - Suction Side')
        plt.plot(x1[541:]/0.0859693,yPlus1[541:],ps,linewidth=2,label = label2 + ' - Pressure Side')

    
    plt.legend(loc=0)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title('Variation of $y^+$ across a T106A blade')
    plt.xlabel('x/C$_{ax}$',fontsize=20)
    plt.ylabel('y$^+$',fontsize=20)
    
def plotWallShearStress(directory,time1,mark,n,label2,fine='no',Re=60000):  
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(n)  
    filename = directory+'surfaces/'+str(time1)+'/wallShearStress_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0] 
    wx  = A[:,3]
    wy = A[:,4]
    wx2  = np.append(wx[0:190], np.append(wx[430:859], np.append(wx[191:429],wx[0])))
    wx1  = np.append(wx2[50:585],np.append(wx2[586:858],wx2[0]))
    wy2  = np.append(wy[0:190], np.append(wy[430:859], np.append(wy[191:429],wy[0])))
    wy1  = np.append(wy2[50:585],np.append(wy2[586:858],wy2[0]))
    
#    x2 = np.append(x[0:190], np.append(x[430:859], np.append(x[191:429],x[0])))
#    x1 = np.append(x2[50:585],np.append(x2[586:858],x2[0]))    
#    plt.plot(x1[:540]/0.0859693,wx1[:540],mark+'--',linewidth=2,label = label2+'wx')
#    plt.plot(x1[:540]/0.0859693,wy1[:540],mark,linewidth=2,label = label2+'wy')
#    
    
    
    
    
    
    norm = 601.56
    if Re==100000:
        norm = 963.5
    wallShearStress = np.sqrt(A[:,3]**2+A[:,4]**2)/norm
    signWSS = np.empty(len(wallShearStress))
    for i in range(len(wallShearStress)):
        if 0<=wx[i]:
            signWSS[i] = -1
        else:
           signWSS[i] = 1     
    wallShearStress = wallShearStress*signWSS
    #this is the normalisation (rho*U^2)_(inlet)
    if fine=='yes':
        wallShearStress2 = np.append(wallShearStress[0:287], np.append(wallShearStress[646:1291], np.append(wallShearStress[288:645],wallShearStress[0]))) #reorder the data around the blade. 
        wallShearStress1 = np.append(wallShearStress2[81:877],np.append(wallShearStress2[878:1291],wallShearStress2[0:82]))
        x2 = np.append(x[0:287], np.append(x[646:1291], np.append(x[288:645],x[0])))
        x1 = np.append(x2[81:877],np.append(x2[878:1291],x2[0:82]))
        plt.plot(x1[:780]/0.0859693,wallShearStress1[:780],mark,linewidth=2,label = label2)
#        plt.plot(x1[811:]/0.859693,yPlus1[811:],ps,linewidth=2,label = label2 + ' - Pressure Side')
    else:
        wallShearStress2 = np.append(wallShearStress[0:190], np.append(wallShearStress[430:859], np.append(wallShearStress[191:429],wallShearStress[0]))) #reorder the data around the blade. 
        wallShearStress1 = np.append(wallShearStress2[50:585],np.append(wallShearStress2[586:858],wallShearStress2[0]))        
        x2 = np.append(x[0:190], np.append(x[430:859], np.append(x[191:429],x[0])))
        x1 = np.append(x2[50:585],np.append(x2[586:858],x2[0]))
        plt.plot(x1[:540]/0.0859693,wallShearStress1[:540],mark,linewidth=2,label = label2)
#        plt.plot(x1[541:]/0.0859693,wallShearStress1[541:],'b',linewidth=2,label = 'Pressure Side')
    
    plt.legend(loc=0)
#    plt.ylim(0,.1)    
    plt.xlabel('x/C$_{ax}$',fontsize=20)
    plt.ylabel('$\\tau_w$',fontsize=20)
    plt.ylim(-0.005,0.04)
    plt.xlim(0,1.01)
    plt.title('Wall Shear Stress')
    plt.grid()
    
def plotCpTimeAv(directory,time1,mark,n,label2,model,fine='no'):   #note that these functions should be combined. 
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(n)  
    filename = directory+'surfaces/'+str(time1)+'/pMean_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0]
    p = A[:,3]
    mixedOutUp = xtractDataTimeAv(directory,time1,'up',model,'MixedOut')
    mixedOutDown = xtractDataTimeAv(directory,time1,'down',model,'MixedOut')
    if fine=='yes':
        p1 = p
        x1 = x
    else:
        p1 = np.append(p[0:190], np.append(p[430:859], p[191:429]))
        x1 = np.append(x[0:190], np.append(x[430:859], x[191:429]))
    cp = (p1- mixedOutDown[0]) / (mixedOutUp[1]-mixedOutDown[0])
    plt.plot(x1/0.859693,cp,mark,linewidth=2,label = label2)
    plt.legend()
    plt.xlim(-0.02,1)
    plt.ylim(-0.7,1.01)
    
    plt.xlabel('x/C$_{ax}$',fontsize=20)
    plt.ylabel('C$_p$',fontsize=20)
    
    
def computeWakeLoss(directory,time1,mark,n,label2,model):   
    import numpy as np
    import matplotlib.pyplot as plt
    mixedOutUp = xtractData(directory,time1,'up',model,'MixedOut')
    mixedOutDown = xtractData(directory,time1,'down',model,'MixedOut')
    Quant = xtractData(directory,time1,'down',model,'wakeLoss')    
    pStag = Quant[2] 
    y = Quant[0]
    wakeLoss = np.empty([1,len(pStag)])
    yStar = -((y - max(y))/(max(y) - min(y)))
    ref = mixedOutUp[1] - mixedOutDown[0]
    wakeLoss  = (mixedOutUp[1] - pStag)/ref
    plt.figure(n)  
    plt.plot(wakeLoss,yStar,mark,linewidth=2,label = label2)
    plt.legend(loc=0)
    plt.xlabel('Omega',fontsize=20)
    plt.ylabel('y$^*$',fontsize=20)
    IntWakeLoss = computeIntegralWakeLoss(wakeLoss,yStar)
    return IntWakeLoss

def computeWakeLossTimeAv(directory,time1,mark,n,label2,model):   
    import numpy as np
    import matplotlib.pyplot as plt
    mixedOutUp = xtractDataTimeAv(directory,time1,'up',model,'MixedOut')
    mixedOutDown = xtractDataTimeAv(directory,time1,'down',model,'MixedOut')
    Quant = xtractDataTimeAv(directory,time1,'down',model,'wakeLoss')    
    pStag = Quant[2] 
    y = Quant[0]
    wakeLoss = np.empty([1,len(pStag)])
    yStar = -((y - max(y))/(max(y) - min(y)))
    ref = mixedOutUp[1] - mixedOutDown[0]
    wakeLoss  = (mixedOutUp[1] - pStag)/ref
    plt.figure(n)  
    plt.plot(wakeLoss,yStar,mark,linewidth=2,label = label2)
    plt.legend(loc=4)
    plt.xlabel('Omega',fontsize=20)
    plt.ylabel('y$^*$',fontsize=20)
    IntWakeLoss = computeIntegralWakeLoss(wakeLoss,yStar)
    return IntWakeLoss

def casesLoading(caseToRun):
    if caseToRun == 'IAVariance':
        k=10
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45/postProcessing/';
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_50/postProcessing/';
        # filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_75/postProcessing/';
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/';
        # filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_25/postProcessing/';
        filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_50/postProcessing/';
        filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/47/postProcessing/';
        
        plotCp(filename1,90,'g--',k,'IA 45$^\circ$','SST')
        plotCp(filename2,90,'c--',k,'IA 45.50$^\circ$','SST')
        # plotCpProfile(filename3,90,'y--',k,'IA 45.75');
        plotCp(filename4,90,'b--',k,'IA 46$^\circ$','SST')
        plotCp(filename5,90,'r--',k,'IA 46.10$^\circ$','SST')
        # plotCpProfile(filename6,90,'c--',k,'IA 46.25');
        plotCp(filename7,90,'y--',k,'IA 46.50$^\circ$','SST')
        plotCp(filename8,90,'m--',k,'IA 47$^\circ$','SST')
        plotTu4LoadingData(k)
        
    elif caseToRun =='TurbModel':
        k=90
#        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelLimitedLinear/postProcessing/'
#        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LimitedLinearSchemes/Tu4/postProcessing/'
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilon/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/realizablekEpsilon/postProcessing/'
#        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilonLimitedLinear/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAUpwind/postProcessing/'
        filename5 ='/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percent/postProcessing/'        
        filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percentStartUpwind/postProcessing/'
        filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percentStartFloor/postProcessing/'
        filename9 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu1percent/postProcessing/'
        filename10 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAQCRLimLinear/postProcessing/'
        filename11 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percentStart/postProcessing/'

#        filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/realizablekEpsilonLimitedLinear/postProcessing/'
        plotCp(filename3,90,'b--',k,'k-$\omega$ SST','SST')
#        plotCp(filename4,60,'r--',k,'Spalart Allmaras','SA')
#        plotCp(filename1,7.994,'g--',k,'k-$\epsilon$','k-Epsilon')
#        plotCp(filename2,88,'c--',k,'Realizable k-$\epsilon$','k-Epsilon')
#        plotCp(filename5,0.3,'m--',k,'Transition 300','SST')
        plotCp(filename5,9.3,'y--',k,'Transition 9300 $','SST')
#        plotCp(filename6,1.69,'g--',k,'Transition 370','SST')
        plotCp(filename5,4.3,'r--',k,'Transition 4300$','SST')
#        plotCp(filename6,3.51,'m--',k,'Transition 3510 LimLinear0.5','SST')
        plotCp(filename8,6.8,'c--',k,'Transition 6800 floor','SST')
        plotCp(filename7,9.009,'b',k,'Transition 3500 Tu1%','SST')
        plotCp(filename10,40,'g',k,'Spalart Allmaras Lim Linear 0.5 QCR','SA')
        plotCp(filename11,4.35,'g',k,'Transition 4350','SST')



        plotTu4LoadingData(k)
        plt.title('Pressure Coefficient on the blade surface at Re$_{2is}$ = 60000')

       
    elif caseToRun == 'Schemes':
        k=78
#        plotCp(filename5,30,'y--',k,'k-$\epsilon$ LimitedLinear','k-Epsilon')
#        plotCp(filename5,54,'y--',k,'k-$\epsilon$ LimitedLinear','k-Epsilon')
#        plotCp(filename7,96,'c--',k,'Realizable k-$\epsilon$ LL','k-Epsilon')
#        plotCp(filename1,60,'r--',k,'SA model LL','SA')
#        plotCp(filename2,30,'b--',k,'SST model LL','SST')
#       #Limited Linear 
    
        filenamea = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilonLimitedLinear/postProcessing/'
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelLimitedLinear/postProcessing/'
        filenameb = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilonLimitedLinear/postProcessing/'
        filenamec = '/home/harshal/OpenFOA9M/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/'
        filename45 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/'
        filenamesul = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/'
        filenameLang = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/'
        lam = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/'
        file0_01 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_01/postProcessing/'
        sst100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/'
        lke60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/'
        lke100k  = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/'
#        plotCp(filename1,60,'r--',k,'SA model LL','SA')
#        plotCp(filenamea,120,'y--',k,'k-9$\epsilon$ LL','k-Epsilon')
#        plotCp(filenameb,200,'c--',k,'Realizable k-$\epsilon$ LL','k-Epsilon')
#        plotCp(filenamec,100,'g--',k,'SST LL','SST')
#        plotCp(filename45,5.792,'m--',k,'SST upwind ','SST')
        plotCp(filename45,74,'r',k,'SST (Turb)','SST')

#        plotCp(filenamesul,50,'b--',k,'Suluksna Trans 50','SST') 
        plotCp(filenamesul,14,'b',k,'$\gamma Re_{\\theta}$ (Trans)','SST') 
        plotCp(filenamesul,14,'b',k,'$\gamma Re_{\\theta}$ (Trans)','SST') 

        plotCp(filenamesul,0.024,'m',k,'URANS $\gamma Re_{\\theta}$ (Trans)','SST') 
#        plotCp(filenamesul,0.024,'m',k,'URANS $\gamma Re_{\\theta}$ (Trans)','SST')
        plotCp(file0_01,40,'brown',k,'URANS $\gamma Re_{\\theta}$ (Trans)_01','SST')
#        plotCp(filenameLang,50,'c--',k,'Lang Trans','SST')  
        plotCp(lam,14,'g',k,'Laminar','SST')         
        plotCp(sst100k,8,'r--',k,'100kSST','SST') 
        plotCp(lke60k,3.2,'m--',k,'60klke 32','SST')  

        plotCp(lke60k,4.8,'r--',k,'60klke 48','SST')  
#        plotCp(lke60k,6,'g--',k,'60klke 6','LKE')  
        plotCp(lke100k,100,'g--',k,'100klke 6.4','LKE')   

        plotTu4LoadingData(k)
        plotTu0LoadingExptData(k) 
#        plotCp(filename3,54,'y--',k,'k-$\epsilon$ LimitedLinear','k-Epsilon')
#        plotCp(filename4,96,'c--',k,'Realizable k-$\epsilon$ LL','k-Epsilon')
        
        
        
    elif caseToRun == 'Laminar':
        k=38
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LaminarCase/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LaminarCaseFeb2107/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinearTry3/postProcessing/'
        plt.figure(k)          
#        plotCp(filename1,30,'g',k,'Laminar 30000','SST')
#        plotCp(filename1,60,'b',k,'Laminar 60000','SST')
#        plotCp(filename1,90,'r',k,'Laminarl 90000','SST')
#        plotCp(filename2,456,'m',k,'Laminar 456k','SST')
        plotCp(filename2,300,'brown',k,'Laminar 300k','SST')
#        plotCp(filename3,0.18,'r',k,'Laminar URANS18','SST')
#        plotCp(filename3,0.34,'grey',k,'Laminar URANS34','SST')
#        plotCp(filename3,0.66,'g',k,'Laminar URANS66','SST')
#        plotCp(filename3,0.74,'b',k,'Laminar URANS74','SST')
        plotCp(filename3,0.9801,'r',k,'Laminar URANS98','SST')
        plotCp(filename3,0.9001,'b',k,'Laminar URANS90','SST')
        plotCp(filename3,0.8201,'c',k,'Laminar URANS82','SST')

#        plotCpTimeAv(filename4,0.022,'c',k,'URANS022','LKE')

        plotTu4LoadingData(k)
        plotTu0LoadingExptData(k)
        
    elif caseToRun == 'URANS':  
        k=46
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Models/SA/SA1/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Models/SA/QCR/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAUpwind/postProcessing/'        
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAQCRLimLinear/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_LESGrid/Steady/Models/SA/SALimLinearQCR/postProcessing/'

        plotCpTimeAv(filename1,0.426,'b--',k,'SA URANS','SA')
        plotCpTimeAv(filename2,0.29,'r--',k,'SA URANS QCR','SA')
#        plotCpTimeAv(filename1,0.456,'c--',k,'T=0.456 URANS','SA')
#        plotCpTimeAv(filename1,0.486,'g--',k,'T=0.486 URANS','SA')
        plotCp(filename4,12.1,'g--',k,'SA RANS Lim Linear 0.5 QCR 12 ','SA')
        plotCp(filename4,20.1,'m--',k,'SA RANS Lim Linear 0.5 QCR 20 ','SA')

#        plotCp(filename3,60,'y--',k,'SA RANS Upwind','SA')
#        plotCp(filename5,150,'c--',k,'SA RANS Lim Linear 0.5 QCR LES grid','SA')

        plotTu4LoadingData(k)
        plotTu0LoadingExptData(k)
        
    elif caseToRun == 'Transition':
        k=48
#        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/30Inlet10/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/35Inlet10/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/40Inlet10/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/43Inlet10/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/46_1Inlet10/postProcessing/'
        filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/'
        filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/55Inlet10/postProcessing/'
        filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/50Inlet10/postProcessing/'
        filenameFine='/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/FinerDNSgrid/TransSteady/Tu4Suluksna/postProcessing/'

##        plotCp(filename1,5,'g--',k,'TM 40Inlet10 ','SST')
#        plotCp(filename2,30,'b--',k,'35$^o$ ','SST')
##        plotCp(filename3,30,'r--',k,'40$^o$ ','SST')
##        plotCp(filename4,30,'g--',k,'43$^o$','SST')
#        plotCp(filename5,12,'c--',k,'46.1$^o$ ','SST')
#        plotCp(filename6,38.8,'y--',k,'LKE 38800','LKE')
#        plotCp(filename6,10.8,'g--',k,'LKE 10800','LKE')
#        plotCp(filename7,12,'r--',k,'55$^o$','SST')
        plotCp(filenameFine,10,'r--',k,'Finer grid $\gamma Re_{\\theta}$','SST','yes')

        plotTu4LoadingData(k)
        plotTu0LoadingExptData(k)   
        
def casesWakeLoss(caseToRun):
    import numpy as np
    import matplotlib.pyplot as plt

    if caseToRun == 'IAVariance':
        k=10
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_50/postProcessing/'
        # filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_75/postProcessing/';
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/'
        # filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_25/postProcessing/';
        filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_50/postProcessing/'
        filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/47/postProcessing/'
        
    elif caseToRun =='TurbModel':
        k=47
        Loss = np.empty([9,1]) #Running here with the upwind scheme. 
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/'       
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilon/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAUpwind/postProcessing/' 
        filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelQCRLimLinear0_5/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/realizablekEpsilon/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percent/postProcessing/'
        filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelQCRUpwind/postProcessing/'
        filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAQCRLimLinear/postProcessing/'
        filename9 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAQCRUpwind/postProcessing/'
        
        
        Loss[0] = computeWakeLoss(filename1,90,'b--',k,'k-$\omega$ SST','SST')
        Loss[1] = computeWakeLoss(filename3,60,'r--',k,'Spalart Allmaras','SA')
        Loss[2] = computeWakeLoss(filename2,7.994,'g--',k,'k-$\epsilon$','k-Epsilon')
        Loss[3] = computeWakeLoss(filename4,88,'c--',k,'Realizable k-$\epsilon$','k-Epsilon')
        Loss[4] = computeWakeLoss(filename5,4.8,'m--',k,'Transition 4800$','SST')
        Loss[5] = computeWakeLoss(filename6,36,'y--',k,'Spalart Allmaras Lim Linear 0.5','SA')
        Loss[6] = computeWakeLoss(filename7,14,'r',k,'Spalart Allmaras Upwind ','SA')
        Loss[7] = computeWakeLoss(filename8,150,'g',k,'Spalart Allmaras Lim Linear 0.5 QCR','SA')
        Loss[8] = computeWakeLoss(filename9,20,'g',k,'Spalart Allmaras Upwind QCR','SA')

#        wake1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/Re60k_data/wake_experiment_clean.dat'
        print Loss
        plotTu4WakeLoss(k)
        plotTu4_Re100kWakeLoss(k)
        plt.legend(loc=4)
#        plotExptWakeLoss(k)
        plt.title('Wake Stagnation Pressure Loss at x=1.26, Re$_{2is}$ = 60000')
        
    elif caseToRun == 'Schemes':
        Loss = np.empty([10,1])
#        Limited Linear scheme - not converged yet. 
        k=47
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelLimitedLinear/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LimitedLinearSchemes/Tu4/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilonLimitedLinear/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/realizablekEpsilonLimitedLinear/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/'
        filename45 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/'
        filenamesul = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/'
        filenameLang = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/'
        lam = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/'
        file0_01 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_01/postProcessing/'
        sst100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/'
        lke100k  = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/'

#        Loss[0] = computeWakeLoss(filename1,100,'r--',k,'SA model LL','SA')
#        Loss[1] = computeWakeLoss(filename4,200,'c--',k,'Realizable k-$\epsilon$','k-Epsilon')
#        Loss[2] = computeWakeLoss(filename3,192.168,'y--',k,'k-$\epsilon$ LL','k-Epsilon')
#        Loss[3] = computeWakeLoss(filename5,100,'g--',k,'SST LL','SST')
#        Loss[0] = computeWakeLoss(filename45,5.792,'m--',k,'SST upwind','SST')
        Loss[1] = computeWakeLoss(filename45,74,'r',k,'SST (Turb)','SST')
        Loss[2] = computeWakeLoss(sst100k,4,'r--',k,'100kSST','SST')  

#        Loss[2] = computeWakeLoss(filenamesul,50,'b--',k,'Suluksna 50 Trans','SST')
        Loss[3] = computeWakeLoss(filenamesul,14,'b',k,'$\gamma Re_{\\theta}$ (Trans)','SST')

#        Loss[4] = computeWakeLoss(filenameLang,50,'c--',k,'Lang Trans','SST')  
        Loss[5] = computeWakeLoss(lam,14,'g',k,'Laminar','SST')     
        Loss[6] = computeWakeLoss(file0_01,40,'brown',k,'URANS $\gamma Re_{\\theta}$ (Trans)_01','SST') 
        Loss[7] = computeWakeLoss(lke100k,100,'g--',k,'100klke 1','LKE')   
        print Loss
        plotTu4WakeLoss(k)
        plotTu4_Re100kWakeLoss(k)
        plt.legend(loc=4)
#        plotExptWakeLoss(k)
        plt.title('Wake Stagnation Pressure Loss at x=1.26$\\times$ C$_{ax}$, Re$_{2is}$ = 60000')
        
    elif caseToRun == 'Laminar':
        k1=380
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LaminarCase/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LaminarCaseFeb2107/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/'

        plt.figure(k1) 
        Loss = np.empty([5,1])
         
#        plotCp(filename1,30,'g',k,'Laminar 30000','SST')
#        plotCp(filename1,60,'b',k,'Laminar 60000','SST')
#        plotCp(filename1,90,'r',k,'Laminar 90000','SST')
        Loss[0] = computeWakeLoss(filename2,456,'m',k1,'Laminar 456k','SST')
        Loss[1] = computeWakeLoss(filename2,300,'brown',k1,'Laminar 300k','SST')
        Loss[2] = computeWakeLoss(filename3,0.9801,'r',k1,'Laminar URANS98','SST')
        Loss[3] = computeWakeLoss(filename3,0.9001,'b',k1,'Laminar URANS90','SST')
        Loss[4] = computeWakeLoss(filename3,0.8201,'c',k1,'Laminar URANS82','SST')

        plotTu4WakeLoss(k1)
        plt.legend(loc=4)
        
    elif caseToRun == 'URANS':  
        k=40
        Loss = np.empty([5,1])
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Models/SA/SA1/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Models/SA/QCR/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAQCRLimLinear/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelLimitedLinear/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_LESGrid/Steady/Models/SA/SALimLinearQCR/postProcessing/'
        Loss[0] = computeWakeLossTimeAv(filename1,0.426,'b--',k,'SA URANS','SA')
        Loss[1] = computeWakeLossTimeAv(filename2,0.29,'g--',k,'SA URANS QCR','SA')
        Loss[2] = computeWakeLoss(filename3,100,'c--',k,'SA RANS','SA')
#        Loss[2] = computeWakeLossTimeAv(filename2,0.25,'m--',k,'SA URANS (Time Averaged QCR2)','SA')

#        Loss[1] = computeWakeLossTimeAv(filename1,0.456,'c--',k,'T=0.456 URANS','SA')
#        Loss[2] = computeWakeLoss(filename4,130,'m--',k,'T=0.486 URANS','SA')

#        Loss[2] = computeWakeLoss(filename4,150,'g--',k,'SA RANS Lim Linear 0.5 QCR','SA')
        Loss[4] = computeWakeLoss(filename4,74,'m--',k,'SA RANS QCR','SA')

#        Loss[4] = computeWakeLoss(filename5,150,'c--',k,'SA RANS Lim Linear 0.5 QCR LES grid','SA')

        print Loss
        plotTu4WakeLoss(k)
        plt.legend(loc=4)
    elif caseToRun == 'Transition':  
        k=50        
        Loss = np.empty([7,1])
        filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/40Inlet10/postProcessing/'
        filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/35Inlet10/postProcessing/'
        filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/46_1Inlet10/postProcessing/'
        filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/46_1Inlet20/postProcessing/'
        filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/46_1Inlet5/postProcessing/'
        filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/'
        filename9 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/'

#        Loss[0] = computeWakeLoss(filename1,30,'g--',k,'TM 40Inlet10 ','SST')
#        Loss[1] = computeWakeLoss(filename2,30,'m--',k,'TM 45Inlet10 ','SST')
        Loss[2] = computeWakeLoss(filename3,12,'b--',k,'TM 46_1Inlet10 ','SST')
#        Loss[3] = computeWakeLoss(filename4,12.5,'b--',k,'TM 46_1Inlet20 ','SST')
#        Loss[4] = computeWakeLoss(filename5,12,'m--',k,'TM 46_1Inlet5 ','SST')
        Loss[5] = computeWakeLoss(filename6,38.8,'r--',k,'LKE38800','LKE')
        Loss[6] = computeWakeLoss(filename6,70,'g--',k,'LKE70000','LKE')
        Loss[3] = computeWakeLoss(filename9,100,'c--',k,'SST LL','SST')

        plotTu4WakeLoss(k)
    
def timeConvergence(directory,model,time):
    """For time convergence we shall check the outlet magnitude of velocity, pressure, k, omega, rho and T (interesting!)"""
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = mpl.cm.winter
    for i in range(0,len(time)):
        field = xtractData(directory,time[i],'down',model,'timeConvergence')
        y = field[0]
        yStar = -((y - max(y))/(max(y) - min(y)))
        label1 = 'Iterations = ' + str(time[i]*1000)
        plt.figure(101)
        plt.plot(field[1],yStar,color=cmap(i / float(len(time))),linewidth=2,label = label1)
        plt.legend()
        plt.xlabel('Pressure [Pa]',fontsize=20)
        plt.ylabel('y$^*$',fontsize=20)
        plt.figure(102)
        plt.plot(field[2],yStar,color=cmap(i / float(len(time))),linewidth=2,label = label1)
        plt.legend()
        plt.xlabel('|Velocity| [m/s]',fontsize=20)
        plt.ylabel('y$^*$',fontsize=20)
        plt.figure(103)
        plt.plot(field[5],yStar,color=cmap(i / float(len(time))),linewidth=2,label = label1)
        plt.legend()
        plt.xlabel('Density [kg/m$^3$]',fontsize=20)
        plt.ylabel('y$^*$',fontsize=20)
        plt.figure(104)
        plt.plot(field[6],yStar,color=cmap(i / float(len(time))),linewidth=2,label = label1)
        plt.legend()
        plt.xlabel('Temperature [K]',fontsize=20)
        plt.ylabel('y$^*$',fontsize=20)




#        plt.figure(2)
#        plt.plot(field[:,0],field[:,2],'r--',linewidth=2,label = 'Velocity')
        
def computeIntegralWakeLoss(wakeLoss, yStar):
    Loss = 0.0    
    for j in range(1,len(yStar)):
        Delta = yStar[j]-yStar[j-1]
        Loss=Loss+((wakeLoss[j]+wakeLoss[j-1])*Delta/2) 
    print Loss
    return Loss
              
def plotTu4LoadingData(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename= '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/DNS_ExptData/Tu_4percent.csv';
    DNS = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(DNS[:,0],DNS[:,1],'k--',linewidth=2,label = 'DNS Tu 4%')
    plt.legend()  
    
def plotTu0LoadingExptData(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/DNS_ExptData/LoadingExptTu_0.csv';
    Expt = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(Expt[:,0],Expt[:,1],'ko',linewidth=2,label = 'Expt Tu 0%')
    plt.legend() 
    
def plotTauwData(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/DNS_ExptData/tauWallTu4.csv';
    Expt = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(Expt[:,0],Expt[:,1],'k--',linewidth=2,label = 'DNS Tu 4%')
    plt.legend() 
    plt.grid()    


def plotTu4WakeLoss(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/DNS_ExptData/WakeLossTu_4_percent.csv';
    Expt = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(Expt[:,0],Expt[:,1],'k--',linewidth=2,label = 'DNS Tu 4%')
    plt.legend(loc=0) 
    computeIntegralWakeLoss(-Expt[:,0],Expt[:,1])
    
def plotExptWakeLoss(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/DNS_ExptData/WakeLossExpt.csv';
    Expt = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(Expt[:,0],Expt[:,1],'ko',linewidth=2,label = 'Exp. (clean)')
    plt.legend() 
    computeIntegralWakeLoss(-Expt[:,0],Expt[:,1])
    
    
def plotTu4_Re100kWakeLoss(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/DNSData/wakeLoss0b4pct.dat';
    A = np.loadtxt(filename)
    plt.figure(k)
    plt.plot(A[:,1],A[:,2],'b--',linewidth=2,label = 'DNS Tu4% Re=100k')
    plt.legend() 
    computeIntegralWakeLoss(A[:,1],A[:,2])