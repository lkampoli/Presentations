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


#for i in range(self.npoints-1):
#            temp1=self.data[i+1,rho_ind]*self.data[i+1,u_ind]
#            temp0=self.data[i,rho_ind]*self.data[i,u_ind]
#            temp_dist=sqrt((self.data[i+1,0]-self.data[i,0])**2
#                          +(self.data[i+1,1]-self.data[i,1])**2
#                          +(self.data[i+1,2]-self.data[i,2])**2)
#        for j in range(5):
#                self.stag_mass_ave[j]+=(temp1*self.stag[i+1,j]+temp0*self.stag[i,j])*temp_dist
#            temp+=(temp1+temp0)*temp_dist
#        for j in range(5):
#            self.stag_mass_ave[j]/=temp
    
#self.stag[i,0]=sqrt(self.data[i,u_ind]**2+self.data[i,v_ind]**2) #vel mag
#            self.stag[i,1]=self.stag[i,0]/sqrt(gamma*self.data[i,p_ind]/self.data[i,rho_ind]) #Mach
#            self.stag[i,2]=self.data[i,T_ind]*(1.+0.2*self.stag[i,1]**2) # stag temp
#            self.stag[i,3]=self.data[i,p_ind]*(1.+0.2*self.stag[i,1]**2)**(1.4/0.4) #stag press
#            self.stag[i,4]=self.data[i,rho_ind]*(1.+0.2*self.stag[i,1]**2)**(1./0.4) #stag density
            
def massAvQuant(rho,p,u,v,y):
    import numpy as np
    from math import sqrt,atan
    gam1=gamma-1.0
    mdot=0.0
    l = np.abs(y[0]-y[-1])
    
    magVel = np.sqrt(u**2+v**2)
    MachNo = magVel/np.sqrt(gamma*p/rho)
#    StagT = self.data[i,T_ind]*(1.+0.2*self.stag[i,1]**2) #can be introduced if required. 
    stagP = p*(1+0.2*MachNo**2)**(gamma/gam1)
    stagRho = rho*(1+0.2*MachNo**2)**(1/gam1)
    
    massAvStagP = 0.0
    massAvP = 0.0 
           
    for j in range(1,len(y)):
        Delta = y[j-1]-y[j]
        mdot=mdot+((rho[j]*u[j]+rho[j-1]*u[j-1])*Delta/2)
        massAvStagP=massAvStagP+ ((rho[j]*u[j]*stagP[j]+rho[j-1]*u[j-1]*stagP[j-1])*Delta/2)
        massAvP=massAvP+ ((rho[j]*u[j]*p[j]+rho[j-1]*u[j-1]*p[j-1])*Delta/2)
    
    
    return[massAvP/mdot,massAvStagP/mdot]

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
            fileOther = directory+'sets/'+str(time1)+'/inlet_p_omega_kt_kl_rho_T.xy'
        elif location =='out':
            fileVelo  = directory+'sets/'+str(time1)+'/outlet_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/outlet_p_omega_kt_kl_rho_T.xy'
        elif location == 'up':
            fileVelo  = directory+'sets/'+str(time1)+'/upstream_U.xy'
            fileOther = directory+'sets/'+str(time1)+'/upstream_p_omega_kt_kl_rho_T.xy'
        elif location == 'down':
            fileVelo   = directory+'sets/'+str(time1)+'/downstream_U.xy'
            fileOther  = directory+'sets/'+str(time1)+'/downstream_p_omega_kt_kl_rho_T.xy'
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
            mixedOutQuant = mixedOutQuantity(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0]) #may need to change this for certain cases
        elif model =='LKE':
            mixedOutQuant = mixedOutQuantity(other[:,5],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        else:
            mixedOutQuant = mixedOutQuantity(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        return mixedOutQuant
    elif purpose == 'MassAvWakeLoss':
        if model == 'SA':
            massAvQuant_1 = massAvQuant(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0]) #may need to change this for certain cases
        elif model =='LKE':
            massAvQuant_1 = massAvQuant(other[:,5],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        else:
            massAvQuant_1 = massAvQuant(other[:,4],other[:,1],velocity[:,1],velocity[:,2],velocity[:,0])
        return massAvQuant_1        
    elif purpose == 'wakeLoss':
        if model == 'SA':
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,4])**0.5
            PStag   = other[:,1]*(1+0.2*MachNo**2)**(1.4/0.4)
        elif model == 'LKE':
            MachNo  = MagVelo/(gamma*other[:,1]/other[:,5])**0.5
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
#    print len(p)
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
    plt.plot(x1/0.0859693,cp,mark, linewidth=2,markersize=5,markevery=20,label = label2)
    plt.legend()
    plt.xlim(-0.02,1)
    plt.ylim(-0.7,1.01)
#    plt.title('Loading on the T106A blade for different inlet angles')
    plt.xlabel('$x/C_{ax}$',fontsize=20)
    plt.ylabel('C$_p$',fontsize=20)
    
def plotyPlus(directory,time1,mark,n,label2,fine='no'):  
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(n)  
    filename = directory+'surfaces/'+str(time1)+'/yPlus_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0]
    yPlus = A[:,3]
#    print yPlus
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
      
    filename = directory+'surfaces/'+str(time1)+'/wallShearStress_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0] 
    wx  = A[:,3]
    wy = A[:,4]
#    
    if fine=='nonDim':
        norm = 1
    elif Re==100000:
        norm = 992.8
    else:
        norm = 596.59 
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
        plt.figure(n)
        plt.plot(x1[:780]/0.0859693,wallShearStress1[:780],mark,marker="o",linewidth=1,markersize=1,markevery=20,label = label2)
        plt.legend(loc=0)
        plt.xlabel('x/C$_{ax}$',fontsize=20)
        plt.ylabel('$\\tau_w$',fontsize=30)
        plt.ylim(-0.005,0.02)
        plt.xlim(0,1.00)
#      plt.title('Wall Shear Stress')
        plt.grid()
        plt.figure(n+2)
        plt.plot(x1[780:]/0.0859693,wallShearStress1[780:],mark,markersize=2,linewidth=2,label = label2)
        plt.legend(loc=0)
        plt.xlabel('x/C$_{ax}$',fontsize=20)
        plt.ylabel('$\\tau_w$',fontsize=30)
        plt.ylim(-0.01,0.04)
        plt.xlim(0,1.0)
        
        
    elif fine == 'nonDim':
        wallShearStress2 = np.append(wallShearStress[0:317], np.append(wallShearStress[718:1435], np.append(wallShearStress[318:717],wallShearStress[0]))) #reorder the data around the blade. 
        wallShearStress1 = np.append(wallShearStress2[84:977],np.append(np.append(wallShearStress2[978:1434],wallShearStress2[0]),wallShearStress2[0:84],)        )
        x2 = np.append(x[0:317], np.append(x[718:1435], np.append(x[318:717],x[0])))
        x1 = np.append(x2[84:977],np.append(np.append(x2[978:1434],x2[0]),x2[0:84]))
        plt.figure(n)
        plt.plot(x[:894]/0.0859693,wallShearStress[:894],mark,linewidth=2,label = label2)
        plt.legend(loc=0)
        plt.xlabel('x/C$_{ax}$',fontsize=20)
        plt.ylabel('$\\tau_w$',fontsize=30)
        plt.ylim(-0.005,0.02)
        plt.xlim(0,1.0)
#        plt.figure(n+2)
#        plt.plot(x1[894:]/0.0859693,wallShearStress1[894:],mark,linewidth=2,label = label2)
#        plt.legend(loc=0)
#        plt.xlabel('x/C$_{ax}$',fontsize=20)
#        plt.ylabel('$\\tau_w$',fontsize=20)
#        plt.ylim(-0.005,0.04)
#        plt.xlim(0,1.01)
#        plt.title('Wall Shear Stress')
        plt.grid()
    
    
         
    else:
        wallShearStress2 = np.append(wallShearStress[0:190], np.append(wallShearStress[430:859], np.append(wallShearStress[191:429],wallShearStress[0]))) #reorder the data around the blade. 
        wallShearStress1 = np.append(wallShearStress2[50:585],np.append(np.append(wallShearStress2[586:858],wallShearStress2[0]),wallShearStress2[0:50],)        )
        x2 = np.append(x[0:190], np.append(x[430:859], np.append(x[191:429],x[0])))
        x1 = np.append(x2[50:585],np.append(np.append(x2[586:858],x2[0]),x2[0:50]))
        plt.figure(n)
        plt.plot(x1[:535]/0.0859693,wallShearStress1[:535],mark, linewidth=2,markersize=5,markevery=20,label = label2)
        plt.legend(loc=0)
        plt.xlabel('x/C$_{ax}$',fontsize=20)
        plt.ylabel('$\\tau_w$',fontsize=30)
        plt.ylim(-0.005,0.04)
        plt.xlim(0,1.0)
        plt.figure(n+2)
        plt.plot(x1[535:]/0.0859693,wallShearStress1[535:],mark,linewidth=2,label = label2)
        plt.legend(loc=0)
        plt.xlabel('x/C$_{ax}$',fontsize=20)
        plt.ylabel('$\\tau_w$',fontsize=20)
        plt.ylim(-0.005,0.04)
        plt.xlim(0,1.0)
#     plt.title('Wall Shear Stress')
        plt.grid()
        
        
def plotWallShearStress_LES(directory,time1,mark,n,label2,Re=60000):  
    import numpy as np
    import matplotlib.pyplot as plt
      
    filename = directory+'surfaces/'+str(time1)+'/wallShearStress_blade.raw'
    A = np.loadtxt(filename)
    x = A[:,0] 
    wx  = A[:,3]
    wy = A[:,4]
#    
    norm = 1
    
    wallShearStress = np.sqrt(A[:,3]**2+A[:,4]**2)/norm
    signWSS = np.empty(len(wallShearStress))
    for i in range(len(wallShearStress)):
        if 0<=wx[i]:
            signWSS[i] = -1
        else:
           signWSS[i] = 1     
    wallShearStress = -wallShearStress*signWSS
    #this is the normalisation (rho*U^2)_(inlet)
    wallShearStress2 = np.append(wallShearStress[49:143], wallShearStress[287:528]) #reorder the data around the blade. 
    wallShearStress1 = np.append(wallShearStress2[84:977],np.append(np.append(wallShearStress2[978:1434],wallShearStress2[0]),wallShearStress2[0:84],)        )
    x2 = np.append(x[49:143], x[287:528]) #reorder the data around the blade. 
    print x2
    x1 = np.append(x2[84:977],np.append(np.append(x2[978:1434],x2[0]),x2[0:84]))
    plt.figure(n)
    plt.plot(x/0.859693,wallShearStress[:],'go',linewidth=2,label = label2)
    plt.legend(loc=0)
    plt.xlabel('x/C$_{ax}$',fontsize=20)
    plt.ylabel('$\\tau_w$',fontsize=30)
    plt.ylim(-0.0005,0.02)
    plt.xlim(0,1.0)
#        plt.figure(n+2)
#        plt.plot(x1[894:]/0.0859693,wallShearStress1[894:],mark,linewidth=2,label = label2)
#        plt.legend(loc=0)
#        plt.xlabel('x/C$_{ax}$',fontsize=20)
#        plt.ylabel('$\\tau_w$',fontsize=20)
#        plt.ylim(-0.005,0.04)
#        plt.xlim(0,1.01)
#        plt.title('Wall Shear Stress')
    plt.grid()
    
    
         
#    else:
#        wallShearStress2 = np.append(wallShearStress[0:190], np.append(wallShearStress[430:859], np.append(wallShearStress[191:429],wallShearStress[0]))) #reorder the data around the blade. 
#        wallShearStress1 = np.append(wallShearStress2[50:585],np.append(np.append(wallShearStress2[586:858],wallShearStress2[0]),wallShearStress2[0:50],)        )
#        x2 = np.append(x[0:190], np.append(x[430:859], np.append(x[191:429],x[0])))
#        x1 = np.append(x2[50:585],np.append(np.append(x2[586:858],x2[0]),x2[0:50]))
#        plt.figure(n)
#        plt.plot(x1[:535]/0.0859693,wallShearStress1[:535],mark,linewidth=2,label = label2)
#        plt.legend(loc=0)
#        plt.xlabel('x/C$_{ax}$',fontsize=20)
#        plt.ylabel('$\\tau_w$',fontsize=30)
#        plt.ylim(-0.005,0.04)
#        plt.xlim(0,1.0)
#        plt.figure(n+2)
#        plt.plot(x1[535:]/0.0859693,wallShearStress1[535:],mark,linewidth=2,label = label2)
#        plt.legend(loc=0)
#        plt.xlabel('x/C$_{ax}$',fontsize=20)
#        plt.ylabel('$\\tau_w$',fontsize=20)
#        plt.ylim(-0.005,0.04)
#        plt.xlim(0,1.0)
##     plt.title('Wall Shear Stress')
#        plt.grid()
          
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
    print ref 
    print mixedOutUp[1]
    print mixedOutDown[0]
    wakeLoss  = (mixedOutUp[1] - pStag)/ref
#    wakeLoss  = pStag/(84.07**2*0.085)

    wakeLoss_MixedOut= (mixedOutUp[1] - mixedOutDown[1])/ref  
#    print 'WakeMixedOutLoss'
#    print label2
#    print wakeLoss_MixedOut
#    print '---'
    plt.figure(n)  
    plt.plot(wakeLoss,yStar,mark, linewidth=2,markersize=5,markevery=10,label = label2)
    plt.legend(loc=0)
    plt.xlabel('$\Omega$',fontsize=20)
#    plt.xlabel('P$_t$',fontsize=20)
    plt.ylabel('$y^*$',fontsize=20)
    IntWakeLoss = computeIntegralWakeLoss(wakeLoss,yStar)
    return [wakeLoss_MixedOut,IntWakeLoss]
    
def computeWakeLoss_massIn(directory,time1,mark,n,label2,model):   
    import numpy as np
    import matplotlib.pyplot as plt
    mixedOutUp = xtractData(directory,time1,'up',model,'MixedOut')
    mixedOutDown = xtractData(directory,time1,'down',model,'MixedOut')
    massAvInlet = xtractData(directory,time1,'up',model,'MassAvWakeLoss')
    massAvOutlet = xtractData(directory,time1,'down',model,'MassAvWakeLoss')
    
#    print massAvInlet[0]
#    (self.inlet.stag_mass_ave[3]-self.outlet.pt_mix_out)/(self.inlet.stag_mass_ave[3]-self.outlet.p_mix_out)
    Quant = xtractData(directory,time1,'down',model,'wakeLoss')    
    pStag = Quant[2] 
    y = Quant[0]
    wakeLoss = np.empty([1,len(pStag)])
    yStar = -((y - max(y))/(max(y) - min(y)))
    ref = massAvInlet[1] - mixedOutDown[0]
#    ref = massAvInlet[0] - massAvOutlet[1]

    wakeLoss  = (massAvInlet[1] - pStag)/ref
    wakeLoss_MassAv  = (massAvInlet[1] - massAvOutlet[1])/ref
#    wakeLoss_MassAv  = (massAvInlet[1] - mixedOutDown[1])/ref

#    print label2    
#    print wakeLoss_MassAv
##    print '***'  
    plt.figure(n)  
    plt.plot(wakeLoss,yStar,mark,linewidth=2,label = label2)
    plt.legend(loc=0)
    plt.xlabel('\Omega',fontsize=20)
    plt.ylabel('y$^*$',fontsize=20)
    IntWakeLoss = computeIntegralWakeLoss(wakeLoss,yStar)
    return [wakeLoss_MassAv,IntWakeLoss]
    

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
    
def plotAll(filename,Timedirec,mark,fig,label,model,Re1=60000,mark2='ro',fine='no'):
#    Loss = np.empty([1,2]) #Running here with the upwind scheme.
    import matplotlib.pyplot as plt
    
    plotCp(filename,Timedirec,mark,fig,label,model)
    Loss = computeWakeLoss(filename,Timedirec,mark,fig+1,label,model)
    
#    print Loss
#    Lossa = computeWakeLoss_massIn(filename,Timedirec,mark,fig+1,label,model)
##    plt.figure(fig+1000)
#    plt.plot(Loss[0],mark2,label=label)
#    plt.plot(Lossa[0],mark2)
#    plt.legend(loc=0)
#    print label
#    print (Loss[0]-Lossa[0])
#    print '-0-0-0-0-0-0-'
    plotWallShearStress(filename,Timedirec,mark,fig+2,label,fine,Re=Re1)
#    print Lossa

def plotDNSData(k):
    import matplotlib as plt
    plotTu4LoadingData(k)
#    plotTu0LoadingExptData(k) 
#    plt.grid()
#    plt.title('Pressure Coefficient on the blade surface at Re$_{2is}$ = 60000')
    plotTu4WakeLoss(k+1)
#    plotExptWakeLoss(k+1)
#    plt.grid()
    plotTauwData(k+2) 
#    plt.grid()
    
def plotDNSData_05(k):
    import matplotlib as plt
    plotTu4LoadingData(k)
#    plotTu0LoadingExptData(k) 
#    plt.grid()
#    plt.title('Pressure Coefficient on the blade surface at Re$_{2is}$ = 60000')
    plotTu4WakeLoss_60ktu05(k+1)
#    plotExptWakeLoss(k+1)
#    plt.grid()
    plotTauwData(k+2) 
#    plt.grid()
    
def plotDNSData100k(k):
    import matplotlib as plt
    plotTu4LoadingData100(k)
#    plotTu0LoadingExptData(k) 
#    plt.grid()
#    plt.title('Pressure Coefficient on the blade surface at Re$_{2is}$ = 60000')
    plotTu4_Re100kWakeLoss(k+1)
#    plt.grid()
    plotTauwDataRe100kTu4(k+2) 
#    plt.grid()
    
def plotLoadingWakeShearStress(caseToRun):
    if caseToRun == 'GammaModelSchemes':
      k=100
      fileTulangtryScheme1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VarySchemes/Tu4LangtryGradOffScheme1/postProcessing/'
      fileTulangtryScheme0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VarySchemes/Tu4LangtryGradOffScheme0/postProcessing/'
      fileTulangtryScheme2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VarySchemes/Tu4LangtryGradOffScheme2/postProcessing/'

      filenameLang = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/'

#      filenamesul = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/'
#      filenameLang = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/'        
        
      plotAll(fileTulangtryScheme1,50,'g',k,'U Scheme linearUpwindV grad(U) ','SST')
      plotAll(fileTulangtryScheme0,50,'b',k,'U Scheme limLinear 0.5','SST')   
      plotAll(fileTulangtryScheme2,50,'m',k,'U Scheme limLinear 1','SST')    

      plotAll(filenameLang,50,'r',k,'U Scheme vanLeerV','SST')   
#      plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/',8,'c--',9, 'SST 100k 8',Re=100000)  
      plotDNSData(k)      
      
    elif caseToRun =='LES1B2U':
        k = 240
        fileLESa = '/media/harshal/GAURA1/IncomingWakes/1B2U/phases/AllPhases/postProcessing/'
#        plotWallShearStress_LES(fileLESa,15,'g',k,'phase15 ','SST')
        fileLESb = '/media/harshal/GAURA1/GEP/LPT_Re60k/TransSteady/kv2Omega/Tu4percentLimLinear/postProcessing/'
        plotWallShearStress_LES(fileLESa,10,'g',k,'phase2 ','SST')
    
    elif caseToRun == 'TransSteadyCompareGamma':
        k = 105
        filenameLang = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/'
        filenamesul = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/'
        

        plotAll(filenameLang,50,'r',k,'Langtry Menter','SST') 
        plotAll(filenamesul,50,'b',k,'Suluksna','SST') 

        plotDNSData(k) 
        
    elif caseToRun == 'LKE60k':
        k =110
        lke60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE60kSteady/Tu4percentLocal/postProcessing/'
        lke60kSteadykklComp2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE60kSteady/Tu4_kklOemgaComp2/postProcessing/'
        
        
        plotAll(lke60k,21.6,'r--',k,'LKE RANS ','LKE') 
#        plotAll(lke60k,12.8,'c--',k,'60klke 12800','LKE') 
        plotAll(lke60k,8,'m--',k,'60klke 8000','LKE') 

        plotAll(lke60k,6.4,'b--',k,'60klke 6400','LKE') 
        plotAll(lke60kSteadykklComp2,20.9,'y--',k,'60k lke 2','LKE') 
 
        plotAll(lke60kSteadykklComp2,23,'c--',k,'60k lke 2','LKE') 
        plotDNSData(k) 
        
    elif caseToRun =='60kTu4':
        k=120
        lam = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/'
        
#       filename45HPTa = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percentHPTa/postProcessing/'

        lke60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE60kSteady/Tu4percentLocal/postProcessing/'
        
        lke60kRANS2_Tu8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE60kSteady/Tu8_kklOemgaComp2/postProcessing/'
        kv2omega = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/kv2Omega/Tu4percentLocal/postProcessing/'
        kv2omegaURANS1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransURANS/kv2Omega/Tu4percentLocal/postProcessing/'
        
        kv2omegaURANS = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/FinalProcess/Links/Tu4URANS_kv2Omega/postProcessing/'
        kv2omegaURANSProd = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransURANS/kv2Omega/Trained/Tu4percentProdTerm/postProcessing/'
        kv2omegaURANSNoProd = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransURANS/kv2Omega/Trained/tu4percentNoProdTerm/postProcessing/'

        kv2omegaNonDim = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/kv2OmegaNonDim/Tu5percent/postProcessing/'
        kv2OmegaLPT60kRect = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases/Re60k/Unsteady/TrainedModelLPtRect60k/postProcessing/'
         
        kv2Tu6perNewProdTerm = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases/Re60k/Unsteady/Tu6percentNewProdTerm/postProcessing/'
        
        k_cont_096 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/OldProdTerm/k_cont_096/postProcessing/'
        
        Rect_kont_096 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WithRamp/Rect_trainedmodel096/postProcessing/'
        Rect_kont_096_prodOff = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WithRamp/RectTrainedModel_ProdBoundOff/postProcessing/'
        kv2OmegaIntAv = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransURANS/kv2Omega/Tu4perLocal_IntAveraging/postProcessing/'
#        plotAll(lam,14,'y',k,'Laminar','SST')      
#        plotAll(filename45HPTa,20,'y',k,'SST - HPT(Turb)','SST')
#        plotAll(kv2omega,25.4,'m',k,'kv2Omega','SST')
#        plotAll(kv2omega,10.8,'m--',k,'kv2Omega 10.8','SST')
#        plotAll(kv2omegaURANS,0.01401,'b',k,'kv2Omega URANS 1401','SST')
#        plotAll(kv2omegaNonDim,6.3,'y--',k,'kv2Omega URANS 6.3','SST',fine='nonDim')
#        plotAll(lke60k,21.6,'r',k,'LKE RANS','LKE') 
#        plotAll(lke60kURANS,0.01401,'y',k,'LKE URANS ','LKE') 
#        plotAll(lke60kRANS2_Tu8,24.3,'b',k,'LKE Tu 8% RANS','LKE')
#        plotAll(filenamesul,50,'b',k,'$\gamma$ Re$_{\\theta}$ RANS','SST') 
        
        filenameSST = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/FinalProcess/Links/Tu4RANS_SST/postProcessing/'
        filenameSA = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SA/Tu4percent/postProcessing/'
        filenameSAQCR = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SAQCR/Tu4percent/postProcessing/'
        filenamesul = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/'
        lke60kURANS = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransURANS/LKE/tu4/postProcessing/'
        kv2NewProdTerm = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/BaseLineCases/Tu4percentNewProdTerm/postProcessing/'
        kv2ModelkCont = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WithRamp/kContour_096_seed0/postProcessing/'
        kwakeOnly = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WakeRegions/kContour_WakeregionOnly/postProcessing/'
        kWake_100kTrained = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WakeRegions/100kTrainedModel_wakeOnly/postProcessing/'
        d60k_100kTrained = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WakeRegions/60k_100kTrained/postProcessing/'
        e60kEnsembleon60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WakeRegions/60kEnsembleon60k/postProcessing/'
        f100kEnsembleon60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WakeRegions/100kEnsembleon60k/postProcessing/'
        g100kRecton60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WithRamp/100kRectOn60k/postProcessing/'
        wallin60k = '/media/harshal/GAURA1/WallinCases/60kWakeRegion/postProcessing/'
#        plotAll(filenameSST,74,'green',k,'k$\omega$-SST','SST')
#        plotAll(filenameSST,74,'m',k,'RANS','SST')

#        plotAll(filenamesul,0.04201,'r',k,'$\gamma$-Re$_{\\theta}$','SST') 
#        plotAll(lke60kURANS,0.02601,'c',k,'LKE 26 ','LKE') 
#        plotAll(lke60kURANS,0.03001,'g',k,'LKE 30 ','LKE') 
###
##        plotAll(kv2omegaURANS,0.01801,'g',k,'kv2${\omega}$ 18','SST')
#        plotAll(kv2omegaURANS,0.02201,'maroon',k,'kv2${\omega}$ 22','SST')
##        plotAll(kv2omegaURANS,0.02601,'b',k,'kv2${\omega}$ 26','SST')
##        plotAll(kv2omegaURANS,0.03001,'pink',k,'kv2${\omega}$ 30','SST')
#        plotAll(kv2omegaURANS,0.03001,'m',k,'kv2${\omega}$ 30','SST')

#        plotAll(kv2omegaURANS,0.03401,'orange',k,'kv2${\omega}$ 34','SST')
#        

#
##        plotAll(kv2omegaURANSProd,0.03502,'b',k,'kv2${\omega}$ Prod Term','SST')
##        plotAll(kv2omegaURANSNoProd,0.04502,'y',k,'kv2${\omega}$ No Prod Term 45','SST')
##        plotAll(kv2omegaURANSNoProd,0.05102,'blue',k,'kv2${\omega}$ No Prod Term 50','SST')
#
##        plotAll(kv2NewProdTerm,0.02501,'orange',k,'kv2${\omega}$ kv2NewProdTerm 25','SST')
##        plotAll(kv2NewProdTerm,0.03001,'purple',k,'kv2${\omega}$ kv2NewProdTerm 30','SST')
##        plotAll(kv2NewProdTerm,0.04501,'springgreen',k,'kv2${\omega}$ kv2NewProdTerm 45','SST')
#        plotAll(kv2NewProdTerm,0.0501,'orange',k,'kv2${\omega}$ kv2NewProdTerm T2','SST')
##        plotAll(kv2Tu6perNewProdTerm,0.03501,'maroon',k,'kv2${\omega}$ 6percent 35','SST')
##        plotAll(kv2Tu6perNewProdTerm,0.04001,'grey',k,'kv2${\omega}$ 6percent 40','SST')
##        plotAll(k_cont_096,0.04901,'pink',k,'k_cont_096','SST')
##        plotAll(Rect_kont_096,0.03501,'orange',k,'kv2${\omega}$ Rect model 35','SST')
##        plotAll(Rect_kont_096,0.04001,'b',k,'kv2${\omega}$ Rect model 40','SST')
##        plotAll(Rect_kont_096,0.04501,'pink',k,'kv2${\omega}$ Rect model 45','SST')
##        plotAll(Rect_kont_096,0.05001,'springgreen',k,'kv2${\omega}$ Rect model 50','SST')
#        plotAll(Rect_kont_096,0.05501,'aqua',k,'kv2${\omega}$ Rect model 55','SST')
#        plotAll(Rect_kont_096,0.06501,'brown',k,'kv2${\omega}$ Rect model 65','SST')
#        plotAll(Rect_kont_096,0.04501,'k',k,'kv2${\omega}$ Rect model prod lim','SST')
#        plotAll(Rect_kont_096,0.04501,'magenta',k,'kv2${\omega}$ Rect model no prod lim','SST')
        
#         plotAll(Rect_kont_096,0.07001,'magenta',k,'kv2${\omega}$ Rect model 70','SST')

#        report figures. 
#        plotAll(filenameSST,74,'g',k,'k${\omega}$-SST','SST')
#        plotAll(filenameSA,100,'orange',k,'SA','SA')
##        plotAll(filenameSAQCR,100,'orange',k,'SA-QCR','SA')
#        plotAll(filenamesul,0.04201,'m',k,'$\gamma$-Re$_{\\theta}$','SST') 
#        plotAll(lke60kURANS,0.03001,'red',k,'LKE','LKE') 
#        plotAll(kv2omegaURANS,0.03801,'blue',k,'k-v$_2$-${\omega}$','SST')
#        plotAll(kv2omegaURANS,0.03801,'blue',k,'Baseline','SST')

#        plotAll(kwakeOnly,0.03401,'m--',k,'k-v$_2$-${\omega}-kwakeonly34$','SST')
#        plotAll(kwakeOnly,0.06001,'g--',k,'k-v$_2$-${\omega}-kwakeonly60$','SST')
        
#        plotAll(Rect_kont_096_prodOff,0.07001,'orange',k,'60k-Rect','SST')
##        plotAll(kv2ModelkCont,0.052501,'m',k,'k-v$_2$-${\omega}$-T2 ','SST')
#        plotAll(kwakeOnly,0.06001,'m',k,'60k-Wake','SST')
##        plotAll(kWake_100kTrained,0.04001,'r',k,'100k-Wake','SST')
#
##        plotAll(kWake_100kTrained,0.04601,'green',k,'100k-Wake','SST')
#        plotAll(kWake_100kTrained,0.07401,'green',k,'100k-Wake','SST')
#        plotAll(d60k_100kTrained,0.05801,'r',k,'60k-100k-Wake','SST')

       
#        plotAll(kv2omegaURANS,0.05801,'maroon',k,'k-v$_2$-${\omega}$ 58','SST')
##        plotAll(kv2OmegaIntAv,0.04401,'pink',k,'k-v$_2$-${\omega}$ 44','SST')
#        
#
#        plotAll(Rect_kont_096,0.07001,'orange',k,'k-v$_2$-${\omega}$-T1','SST')
##        plotAll(kv2ModelkCont,0.024501,'r',k,'kv2${\omega}$ kcont 24','SST')
#
#        plotAll(kv2ModelkCont,0.026501,'m',k,'kv2${\omega}$ Kcont 26','SST')
#        plotAll(kv2ModelkCont,0.034501,'y',k,'kv2${\omega}$ Kcont 34','SST')


##        plotAll(Rect_kont_096_prodOff,0.06501,'r',k,'kv2${\omega}$ Trained ProdBoundoff 65','SST')
#
#        plotAll(kv2NewProdTerm,0.0501,'blue',k,'k-v$_2$-${\omega}$','SST')



# Paper figures - 1: URANS of models 

#        plotAll(filenameSST,74,'go-',k,'k${\omega}$-SST','SST')
#        plotAll(filenameSA,100,'orange',k,'SA','SA')
##        plotAll(filenamesul,0.04201,'ms-',k,'$\gamma$-Re$_{\\theta}$','SST') 
###        plotAll(filenamesul,0.04201,'g',k,'RANS','SST') 
#        plotAll(lke60kURANS,0.02601,'gx-',k,'LKE','LKE') 
#
##        plotAll(lke60kURANS,0.03001,'rx-',k,'LKE','LKE') 
#        plotAll(kv2omegaURANS,0.03801,'b^-',k,'k-v$_2$-${\omega}$','SST') 


#Paper figures - 2: Trained   
        plotAll(kv2omegaURANS,0.03801,'b^-',k,'Baseline','SST')
#        plotAll(Rect_kont_096_prodOff,0.07001,'orange',k,'60k-Rect','SST')
#        plotAll(g100kRecton60k,0.06201,'c*-',k,'100k-Rect','SST')
        plotAll(kwakeOnly,0.06001,'ms-',k,'60k-Wake','SST')
#        plotAll(kWake_100kTrained,0.07401,'go-',k,'100k-Wake','SST')
#        plotAll(d60k_100kTrained,0.05801,'rx-',k,'60k-100k-Wake','SST')
#        plotAll(e60kEnsembleon60k,0.05801,'gs-',k,'60k-Wake-E','SST') #60k-Wake-E
#        plotAll(f1010kEnsembleon60k,0.06001,'ro-',k,'100k-Wake-E','SST')
        plotAll(wallin60k,0.06661,'yp-',k,'Wallin & Johansson','SST')

        



        plotDNSData(k) 
        
    elif caseToRun =='60kTu05':
        k=1233
        kv2_tu075actual = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransURANS/kv2Omega/Tu075percentActual/postProcessing/'
        TrainedRect_60ktu4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/60k/UnSteady/WithRamp/Tu075percent/postProcessing/'       
#        plotAll(kv2_tu075actual,0.04601,'blue',k,'k-v$_2$-${\omega}$','SST')
        plotAll(kv2_tu075actual,0.06401,'b',k,'k-v$_2$-${\omega}$','SST')

#        plotAll(TrainedRect_60ktu4,0.04801,'red',k,'k-v$_2$-${\omega}-Rect-60ktu4$','SST')
#        plotAll(TrainedRect_60ktu4,0.06601,'green',k,'k-v$_2$-${\omega}-Rect-60ktu4$','SST')
        plotAll(TrainedRect_60ktu4,0.08401,'orange',k,'k-v$_2$-${\omega}-Rect-60ktu4$','SST')
       
        plotDNSData_05(k) 
        
        
        
    elif caseToRun == '100kTu4':
        k=1350
        sst100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/'
        lke100k  = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/'
        lke100k_01URANS = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransURANS/tu4/postProcessing/'
        gamma100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/Tu4LangtryGradOff/postProcessing/'
        SA100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SA/postProcessing/'
        lke100ktu6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu6percent/postProcessing/'
        lke100tu4_kkl2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu6percent/postProcessing/'
        kv2OmegaNewProd = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/BaseLineCases/Tu4percentNewProdTerm100k/postProcessing/'
        kv2Omega100k_kcontour_seed0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/100k_kcontour_seed0/postProcessing/'
        
        kv2OmegaOrigProdTerm = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/BaseLineCases/Tu4percent100kOrigProdTerm/postProcessing/'        
        kv2Omega60kRect =  '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/60kTrainedModel_096/postProcessing/'
        kv2OmegakCont = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/60k_kcontour/postProcessing/'
        
        a100kWakeOnly = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/WakeOnly/Trained100kWake/postProcessing/'
        b100k_with60ktrainedWakeonly = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/WakeOnly/Trained60kWake/postProcessing/'
        c60k_100k_trainedWake = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/WakeOnly/60k_100kTrained/postProcessing/'
        d100kEnsembleOn100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/WakeOnly/100kEnsembleon100k/postProcessing/'
        e60kEnsembleOn100k =  '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/WakeOnly/60kEnsembleOn100k/postProcessing/'
        f100kRecton100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/kv2OmegaCases_009/100k/Unsteady/100kRectOn100k/postProcessing/'
        wallin100k = '/media/harshal/GAURA1/WallinCases/100kWakeRegion/postProcessing/'
#        plotAll(sst100k,8,'g',k,'SST (Turb)','SST',100000) 
#        plotAll(lke100k,100,'r',k,'LKE RANS','LKE',100000)   
#        plotAll(lke100ktu6,8.4,'y',k,'LKE RANS Tu6','SST',100000)   
#        plotAll(lke100k_01URANS,0.02601,'m',k,'LKE URANS','SST',100000)
        
#        plotAll(lke100tu4_kkl2,8.4,'c',k,'LKE URANS kkl2','LKE',100000)
#        plotAll(kv2OmegaNewProd,0.02001,'orange',k,'LKE kv2omega newProd 20','SST',100000)
#        plotAll(kv2OmegaNewProd,0.03501,'gold',k,'LKE kv2omega newProd 35','SST',100000)
#        plotAll(kv2OmegaNewProd,0.04001,'purple',k,'LKE kv2omega newProd 40','SST',100000)
#        plotAll(kv2OmegaNewProd,0.05501,'aqua',k,'LKE kv2omega newProd 55','SST',100000)
#        plotAll(kv2OmegaNewProd,0.0601,'springgreen',k,'LKE kv2omega newProd 60','SST',100000)
#        plotAll(kv2Omega60kRect,0.03001,'b',k,'kv2100k_60ktrained 30','SST',100000)
#        plotAll(kv2Omega60kRect,0.06501,'pink',k,'kv2100k_60ktrained 65','SST',100000)
#        plotAll(kv2OmegaOrigProdTerm,0.06401,'b',k,'k-v$_2$-${\omega}$','SST',100000)
#        plotAll(kv2OmegaOrigProdTerm,0.06401,'b',k,'Baseline','SST',100000)
#
#
##        plotAll(kv2Omega60kRect,0.07001,'orange',k,'k-v$_2$-$\omega$-60T1','SST',100000)
##        plotAll(kv2OmegakCont,0.04501,'m',k,'k-v$_2$-$\omega$-60T2','SST',100000)
##        plotAll(kv2Omega100k_kcontour_seed0,0.03601,'green',k,'k-v$_2$-$\omega$-100T2-36','SST',100000)
##        plotAll(kv2Omega100k_kcontour_seed0,0.06201,'green',k,'k-v$_2$-$\omega$-100T2','SST',100000)
##        plotAll(a100kWakeOnly,0.03201,'r--',k,'k-v$_2$-$\omega$-100-100kwakeonly','SST',100000)
##        plotAll(a100kWakeOnly,0.05001,'green',k,'k-v$_2$-$\omega$-100-100kwakeonly','SST',100000)
##        plotAll(a100kWakeOnly,0.07601,'maroon',k,'k-v$_2$-$\omega$-100-60kwakeonly76','SST',100000)
#
##        plotAll(b100k_with60ktrainedWakeonly,0.03001,'pink',k,'k-v$_2$-$\omega$-100-60kwakeonly','SST',100000)
##        plotAll(b100k_with60ktrainedWakeonly,0.05001,'purple',k,'k-v$_2$-$\omega$-100-60kwakeonly50','SST',100000)



       
#        
##        plotAll(kv2OmegaOrigProdTerm,0.03401,'c',k,'kv2OmegaOrig 34','SST',100000)
#        plotAll(kv2OmegaOrigProdTerm,0.04001,'r',k,'kv2OmegaOrig 40','SST',100000)

        #100k report 1
#        plotAll(sst100k,8,'go-',k,'k$\omega$-SST','SST',100000) 
#        plotAll(SA100k,100,'orange',k,'SA','SA',100000)            
#        plotAll(gamma100k,8,'ms-',k,'$\gamma$-Re$_{\\theta}$','SST',100000)
#        plotAll(lke100k_01URANS,0.02601,'rx-',k,'LKE','SST',100000)
#        plotAll(kv2OmegaOrigProdTerm,0.04001,'b^-',k,'k-v$_2$-${\omega}$','SST',100000)
#        plotAll(kv2OmegaOrigProdTerm,0.08201,'y',k,'k-v$_2$-${\omega} 82$','SST',100000)

        #100kReport2
        plotAll(kv2OmegaOrigProdTerm,0.06401,'b^-',k,'Baseline','SST',100000)
#        plotAll(kv2Omega60kRect,0.07001,'orange',k,'60k-Rect','SST',100000)
#        plotAll(f100kRecton100k,0.06001,'c*-',k,'100k-Rect','SST',100000)
#        plotAll(b100k_with60ktrainedWakeonly,0.05001,'ms-',k,'60k-Wake','SST',100000)
        plotAll(a100kWakeOnly,0.05001,'go-',k,'100k-Wake','SST',100000)
#        plotAll(c60k_100k_trainedWake,0.05401,'rx-',k,'60k-100k-Wake','SST',100000) 
#        plotAll(e60kEnsembleOn100k,0.06401,'gs-',k,'60k-Wake-E','SST',100000)         
#        plotAll(d100kEnsembleOn100k,0.06401,'ro-',k,'100k-Wake-E','SST',100000) 

        plotAll(wallin100k,0.1301,'yp-',k,'Wallin & Johansson','SST',100000) 

        plotDNSData100k(k) 


    elif caseToRun == 'LKE100k':
        k =115
        lke60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/'
        lke100k  = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/'
        lke100k_005 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/tu4_005/postProcessing/'
        lke100k_001 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/tu4_001/postProcessing/'
        lke100k_005_tu0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/tu0_1_005/postProcessing/'
        sst100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/'
        lke100k_01URANS = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransURANS/tu4/postProcessing/'
#        plotAll(lke60k,21.6,'r--',k,'60klke 21600','LKE') 
#        plotAll(lke60k,12.8,'c--',k,'60klke 12800','LKE') 
#        plotAll(lke60k,8,'m--',k,'60klke 8000','LKE') 

#        plotAll(lke60k,6.4,'b--',k,'60klke 6400','LKE') 
        plotAll(lke100k,100,'r',k,'LKE RANS L = 10% C$_{ax} Re=100k$','LKE',100000)   
#        plotAll(lke100k_005,100,'r--',k,'100klke L = 5% C$_{ax}$','LKE',100000)   
#        plotAll(lke100k_001,120,'c--',k,'100klke L = 1% C$_{ax}$','LKE',100000)   
#        plotAll(lke100k_005_tu0,100,'m--',k,'100klke L = 5% C$_{ax} Tu=0%$','LKE',100000)
#        plotAll(lke100k_01URANS,0.01801,'g--',k,'100klke URANS L = 10% C$_{ax} 1801 $','LKE',100000)   
        plotAll(lke100k_01URANS,0.02601,'b',k,'LKE URANS L = 10% C$_{ax} Re=100k$','LKE',100000)   
#        plotAll(lke100k_01URANS,0.03401,'c--',k,'LKE URANS L = 10% C$_{ax} 3401$','LKE',100000)   


        plotAll(sst100k,8,'green',k,'k$\omega$-SST','SST',100000) 



        plotDNSData100k(k)         
    elif caseToRun =='TurbModel':
        k=90 # SA, SA=QCR, SST
        filename45 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/'
        lam = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/'


        plotAll(filename45,74,'g',k,'SST (Turb)','SST')
        plotAll(lam,14,'b',k,'Laminar','SST')         
        plotDNSData(k)      
       
    elif caseToRun == 'Schemes':
        k=78
#        
    
        filename45 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/'
        filenamesul = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/'
        filenameLang = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/'
        lam = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/'
        file0_01 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_01/postProcessing/'
        sst100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/'
        lke60k = '/hom60kTu4e/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/'
        lke100k  = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/'
        
        
        
        
        Loss[1] = computeWakeLoss(filename45,74,'r',k,'SST (Turb)','SST')
        Loss[2] = computeWakeLoss(sst100k,4,'r--',k,'100kSST','SST')  

#        Loss[2] = computeWakeLoss(filenamesul,50,'b--',k,'Suluksna 50 Trans','SST')
        Loss[3] = computeWakeLoss(filenamesul,14,'b',k,'$\gamma Re_{\\theta}$ (Trans)','SST')

#        Loss[4] = computeWakeLoss(filenameLang,50,'c--',k,'Lang Trans','SST')  
        Loss[5] = computeWakeLoss(lam,14,'g',k,'Laminar','SST')     
        Loss[6] = computeWakeLoss(file0_01,40,'brown',k,'URANS $\gamma Re_{\\theta}$ (Trans)_01','SST') 
        Loss[7] = computeWakeLoss(lke100k,100,'g--',k,'100klke 1','LKE')   
        Loss[8] = computeWakeLoss(lke60k,6.4,'r--',k,'60klke 48','LKE')  
#        plotCp(filename1,60,'r--',k,'SA model LL','SA')
#        plotCp(filenamea,120,'y--',k,'k-9$\epsilon$ LL','k-Epsilon')
#        plotCp(filenameb,200,'c--',k,'Realizable k-$\epsilon$ LL','k-Epsilon')
#        plotCp(filenamec,100,'g--',k,'SST LL','SST')
#        plotCp(filename45,5.792,'m--',k,'SST upwind ','SST')
        plotCp(filename45,74,'r',k,'SST (Turb)','SST')

#        plotCp(filenamesul,50,'b--',k,'Suluksna Trans 50','SST') 
        plotCp(filenamesul,14,'b',k,'$\gamma Re_{\\theta}$ (Trans)','SST') 

        plotCp(filenamesul,0.024,'m',k,'URANS $\gamma Re_{\\theta}$ (Trans)','SST') 
#        plotCp(filenamesulplotDNSData_05,0.024,'m',k,'URANS $\gamma Re_{\\theta}$ (Trans)','SST')
        plotCp(file0_01,40,'brown',k,'URANS $\gamma Re_{\\theta}$ (Trans)_01','SST')
#        plotCp(filenameLang,50,'c--',k,'Lang Trans','SST')  
        plotCp(lam,14,'g',k,'Laminar','SST')         
        plotCp(sst100k,8,'r--',k,'100kSST','SST') 
        plotCp(lke60k,3.2,'m--',k,'60klke 32','SST')  

#        plotCp(lke60k,12.8,'r--',k,'60klke 48','LKE')  
#        plotCp(lke60k,6,'g--',k,'60klke 6','LKE')  
        plotCp(lke100k,100,'g--',k,'100klke 6.4','LKE')   

        plotTu4LoadingData(k)
#        plotTu0LoadingExptData(k) 
#        plotCp(filename3,54,'y--',k,'k-$\epsilon$ LimitedLinear','k-Epsilon')
#        plotCp(filename4,96,'c--',k,'Realizable k-$\epsilon$ LL','k-Epsilon')
        
        
        
    elif caseToRun == '100k':
        k=38
        sst100k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/'
        lke100k  = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/'
#        plotCpTimeAv(filename4,0.022,'c',k,'URANS022','LKE')

        plotTu4LoadingData(k)
        plotTu0LoadingExptData(k)
        
    elif caseToRun == 'URANS':  
        k=46
        

        plotTu4LoadingData(k)
        plotTu0LoadingExptData(k)
        
    
        
def casesWakeLoss(caseToRun):
    import numpy as np
    import matplotlib.pyplot as plt

    if caseToRun == 'LangM':
        k=10
        
        
    elif caseToRun =='TurbModel':
        k=47
        Loss = np.empty([9,1]) #Running here with the upwind scheme. 
        print Loss
        
        
        
#       
        plotTu4WakeLoss(k)
        plotTu4_Re100kWakeLoss(k)
        plt.legend(loc=4)
#        plotExptWakeLoss(k)
#        plt.title('Wake Stagnation Pressure Loss at x=1.26, Re$_{2is}$ = 60000')
        
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
        lke60k = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/'

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
        Loss[8] = computeWakeLoss(lke60k,6.4,'r--',k,'60klke 48','LKE')  

        print Loss
        plotTu4WakeLoss(k)
        plotTu4_Re100kWakeLoss(k)
        plt.legend(loc=4)
#        plotExptWakeLoss(k)
        plt.title('Wake Stagnation Pressure Loss at x=1.26$\\times$ C$_{ax}$, Re$_{2is}$ = 60000')
        
    elif caseToRun == 'Laminar':
        
       
        plotTu4WakeLoss(k1)
        plt.legend(loc=4)
        
    elif caseToRun == 'URANS':  
        k=40
        

        print Loss
        plotTu4WakeLoss(k)
        plt.legend(loc=4)
    elif caseToRun == 'Transition':  
        k=50        
        Loss = np.empty([7,1])
        

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
    import numpy as np
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/60k/0b_4pct_time_mean/blade.dat';
    cd = np.loadtxt(filename,dtype='float',skiprows=1)
    plt.figure(k)
    plt.plot(cd[:,0]/0.859693,cd[:,4],'k--',linewidth=2,label = 'DNS')
#    plt.legend(loc=0) 
    plt.legend(loc=1,fontsize=15) 
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
#    plt.xlabel(fontsize=20)
#    plt.ylabel(fontsize=20)
    
    plt.title('Pressure Coefficient',fontsize=15)
#    plt.title('Re$_{2is}$=60000, L$_{in}$ = 5% C$_{ax}$, Tu$_{in}$ = 4%')
    plt.grid()
#    plt.savefig(dirPath+'Wake/magVel/magVel_'+str(titleWake[i])+'.jpeg',dpi=500)
    
def plotTu4LoadingData100(k):
    import matplotlib.pyplot as plt
    import numpy as np
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/100k/0b_4pct_time_mean/blade.dat';
    cd = np.loadtxt(filename,dtype='float',skiprows=1)
    plt.figure(k)
    plt.plot(cd[:,0]/0.859693,cd[:,4],'k--',linewidth=2,label = 'DNS')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc=0) 
#    plt.title('Re$_{2is}$=100000, L$_{in}$ = 5% C$_{ax}$, Tu$_{in}$ = 4%')
    plt.grid()   
    
     
def plotTu0LoadingExptData(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/DNS_ExptData/LoadingExptTu_0.csv';
    Expt = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(Expt[:,0],Expt[:,1],'ko',linewidth=2,label = 'Expt Tu 0%')
    plt.legend(loc=0) 
#    plt.title('Re$_{2is}$=60000, L = 10% C$_{ax}$, Tu = 4%')
    
def plotTauwData(k):
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/60k/0b_4pct_time_mean/blade.dat';
    cdr = np.loadtxt(filename,dtype='float',skiprows=1)
    A = cdr[:,0]
    B = cdr[:,-1]
    plt.figure(k)

    plt.plot(cdr[325:,0]/0.859693,cdr[325:,-1],'k--',linewidth=2,label = 'DNS')

#    plt.plot(tauSS[:,0],tauSS[:,1],'k--',linewidth=2,label = 'DNS')

    file2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/60k/0b_4pct_time_mean/tauwall60kSS.csv'
    
    
    tauSS = genfromtxt(file2,delimiter=',')
    plt.figure(k)
#    plt.plot(tauSS[:,0],tauSS[:,1],'k--',linewidth=2,label = 'DNS')
    plt.legend(loc=0,fontsize=15) 
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('Wall Shear Stress - Suction Side',fontsize=15)
#    plt.xlabel(fontsize=20)
#    plt.ylabel(fontsize=20)
    plt.xlim(0,1)
#    plt.legend(loc=0) 
#    plt.title('Re$_{2is}$=60000, L = 5% C$_{ax}$, Tu = 4%')
    plt.grid()   
    
    plt.figure(k+2)
    plt.plot(A[:325]/0.859693,-B[:325],'k--',linewidth=2,label = 'DNS')
    plt.legend(loc=0) 
#    plt.title('Re$_{2is}$=60000, L = 5% C$_{ax}$, Tu = 4%')
    plt.grid()  
    
#    plt.plot(x1[535:850]/0.0859693,wallShearStress1[535:850],mark,linewidth=2,label = label2)
    
def plotTauwDataRe100kTu4(k):
    import matplotlib.pyplot as plt
    import numpy as np
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/100k/0b_4pct_time_mean/blade.dat';
    cd = np.loadtxt(filename,dtype='float',skiprows=1)
    plt.figure(k)
    plt.plot(cd[500:,0]/0.859693,cd[500:,-1],'k--',linewidth=2,label = 'DNS')
#    plt.title('Re$_{2is}$=100000, L$_{in}$ = 5% C$_{ax}$, Tu$_{in}$ = 4%')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plotTauwData(k)
    plt.legend(loc=0) 

    plt.grid() 
    
def plotTauwDataRe100kTu0(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/100ktu0ShearStress.csv';
    Expt = genfromtxt(filename,delimiter=',')
    plt.figure(k)
    plt.plot(Expt[:,0],Expt[:,1],'k--',linewidth=2,label = 'DNS')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc=0) 
#    plt.title('Re$_{2is}$=100000, L$_{in}$ = 5% C$_{ax}$, Tu$_{in}$ = 4%')
    plt.grid() 


def plotTu4WakeLoss_60ktu05(k):
    import matplotlib.pyplot as plt
    import numpy as np
    filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/DNS_ExptData/dns0_5_60k.csv'; # tu0ExptLoss60k_try2.
    Expt2 = genfromtxt(filename2,delimiter=',')
    plt.figure(k)
    plt.plot(Expt2[:,0],-Expt2[:,1]+1,'k--',linewidth=2,label = 'DNS')
    plt.legend(loc=0) 


def plotTu4WakeLoss(k):
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import genfromtxt

    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/60k/0b_4pct_time_mean/outlet.dat';
    cd = np.loadtxt(filename,dtype='float',skiprows=1)
    plt.figure(k)
    plt.plot(cd[:,7],cd[:,1],'k--',linewidth=2,label = 'DNS')

    computeIntegralWakeLoss(-cd[:,7],cd[:,1])    

    
    filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/DNS_ExptData/dns0_5_60k.csv'; # tu0ExptLoss60k_try2.
    filpt_DNS40 = '/home/harshal/PLATUS/harshal/pstag40.csv'
    pt_0 = genfromtxt(filpt_DNS40,delimiter=',')
    pt = np.array([pt_0[1:,0], pt_0[1:,-1]]).transpose()
#    filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/60k/0b_4pct_time_mean/Michel2014_tu4_fred0.csv'
    Expt2 = genfromtxt(filename2,delimiter=',')
    plt.figure(k)
#    plt.plot(pt[:,0],pt[:,1],'k--',linewidth=2,label = 'DNS')
    plt.legend(loc=0) 
#    computeIntegralWakeLoss(-Expt2[:,0],Expt2[:,1])    
#    
#    a100kDNS40_0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/Re100k_data/60k_40chordDown_DNS.csv'
#    a100kDNS30_0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/Re100k_data/60k_30chordUp_DNS.csv'
#    a100kDNS40_down = genfromtxt(a100kDNS40_0,delimiter=',')  
#    a100kDNS30_up = genfromtxt(a100kDNS30_0,delimiter=',')  
#    
#    a100kDNS40_down1 = np.array(a100kDNS40_down[1:,:])#*0.09*0.09
#    a100kDNS30_up1 = np.array(a100kDNS30_up[1:,:])#*0.09*0.09
#    
#    #Generate 100k DNS loss profile. 
#    mixedDown = mixedOutQuantity(a100kDNS40_down1[:,1],a100kDNS40_down1[:,5],a100kDNS40_down1[:,2],a100kDNS40_down1[:,3],a100kDNS40_down1[:,0])
#    mixedUp = mixedOutQuantity(a100kDNS30_up1[:,1],a100kDNS30_up1[:,5],a100kDNS30_up1[:,2],a100kDNS30_up1[:,3],a100kDNS30_up1[:,0])
#    MagVeloDown = np.sqrt(a100kDNS40_down1[:,2]**2+a100kDNS40_down1[:,3]**2)
#    MachNo = MagVeloDown/(1.4*a100kDNS40_down1[:,5]/a100kDNS40_down1[:,1])**0.5
#    PStag   = a100kDNS40_down1[:,5]*(1+0.2*MachNo**2)**(1.4/0.4)
#    y = a100kDNS40_down1[:,0]
#    yStar = -((y - max(y))/(max(y) - min(y)))
#    ref = mixedUp[1] - mixedDown[0]
#    wakeLoss  = (mixedUp[1] - PStag)/ref
#    computeIntegralWakeLoss(-wakeLoss,yStar)
#    plt.figure(k)
#    plt.plot(wakeLoss,yStar,'k--',linewidth=2,label = 'DNS')    
    
    
    
    plt.legend(loc=0,fontsize=15) 
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
#    plt.title('Wake Stag. Pressure Loss, 40% C Downstream, Tu$_{in}$ = 0.5%',fontsize=15)
    plt.xlabel('$\Omega$')
#    plt.title('Re$_{2is}$=60000, L$_{in}$ = 5% C$_{ax}$, Tu$_{in}$ = 4%')
    plt.grid()     
    
   
def plotExptWakeLoss(k):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
#    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/DNS_ExptData/WakeLossExpt.csv';
#    Expt = genfromtxt(filename,delimiter=',')
#    plt.figure(k)
#    plt.plot(Expt[:,0],Expt[:,1],'ko',linewidth=2,label = 'Exp. (clean)')    
    filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/DNS_ExptData/expt0_5_60k.csv'; # tu0ExptLoss60k_try2.
    Expt2 = genfromtxt(filename2,delimiter=',')
    plt.figure(k)
    plt.plot(Expt2[:,0],Expt2[:,1],'go',ms=8,label = 'Expt.')
    plt.legend(loc=0) 
    computeIntegralWakeLoss(-Expt2[:,0],Expt2[:,1])
    
    
def plotTu4_Re100kWakeLoss(k):
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import genfromtxt

    filename = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/100k/0b_4pct_time_mean/outlet.dat';
    cd = np.loadtxt(filename,dtype='float',skiprows=1)
    cd1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/PostProcessData/100k/100kLoss.csv'
    cd0 = genfromtxt(cd1,delimiter=',')
#    print cd0[:,0]
    plt.figure(k)
    plt.plot(cd0[:,1],cd0[:,0],'k--',linewidth=2,label = 'DNS')
    print 'DNS100k Integral WakeLoss'
    computeIntegralWakeLoss(-cd0[:,1],cd0[:,0])

#    a100kDNS40_0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/Re100k_data/100k_40chordDown_DNSmod.csv'
#    a100kDNS30_0 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/PostProcessData/Re100k_data/100k_30chordUp_DNS_mod.csv'
#    a100kDNS40_down = genfromtxt(a100kDNS40_0,delimiter=',')  
#    a100kDNS30_up = genfromtxt(a100kDNS30_0,delimiter=',')  
#    
#    a100kDNS40_down1 = np.array(a100kDNS40_down[1:,:])#*0.09*0.09
#    a100kDNS30_up1 = np.array(a100kDNS30_up[1:,:])#*0.09*0.09
#    
#    #Generate 100k DNS loss profile. 
#    mixedDown = mixedOutQuantity(a100kDNS40_down1[:,1],a100kDNS40_down1[:,5],a100kDNS40_down1[:,2],a100kDNS40_down1[:,3],a100kDNS40_down1[:,0])
#    mixedUp = mixedOutQuantity(a100kDNS30_up1[:,1],a100kDNS30_up1[:,5],a100kDNS30_up1[:,2],a100kDNS30_up1[:,3],a100kDNS30_up1[:,0])
#    MagVeloDown = np.sqrt(a100kDNS40_down1[:,2]**2+a100kDNS40_down1[:,3]**2)
#    MachNo = MagVeloDown/(1.4*a100kDNS40_down1[:,5]/a100kDNS40_down1[:,1])**0.5
#    PStag   = a100kDNS40_down1[:,5]*(1+0.2*MachNo**2)**(1.4/0.4)
#    y = a100kDNS40_down1[:,0]
#    yStar = -((y - max(y))/(max(y) - min(y)))
#    ref = mixedUp[1] - mixedDown[0]
#    wakeLoss  = (mixedUp[1] - PStag)/ref
#    plt.figure(k)
#    plt.plot(wakeLoss,yStar,'k--',linewidth=2,label = 'DNS')
        
    plt.xlabel('$\Omega$')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)       
#    plt.plot(cd[:,7],cd[:,1],'k--',linewidth=2,label = 'DNS')
    plt.legend(loc=0,fontsize=15) 
#    plt.title('Re$_{2is}$=100000, L$_{in}$ = 5% C$_{ax}$, Tu$_{in}$ = 4%')
    plt.grid()    