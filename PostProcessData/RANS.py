# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:46:50 2016
@author: harshal
"""
import PostProcess
from numpy import genfromtxt
import matplotlib.pyplot as plt
import numpy as np
#plt.close("all")
Loading = 'Schemes'
wakeLoss = 'Schemes'
runCases = '60kTu4'  
#runCases = 'LES1B2U'
time = [30,60,90] 
#casesLoading(Loading)
#casesWakeLoss(wakeLoss)

PostProcess.plotLoadingWakeShearStress(runCases)

#plotAll(filenameSST,74,'green',k,'k$\omega$-SST','SST')
#filenameSST = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/FinalProcess/Links/Tu4RANS_SST/postProcessing/'
#plotyPlus(filenameSST,74,'g',908,'k$\omega$-SST','SST')

#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',52.4,'m',5,'LKE 52400')
#plotyPlus('/home/harshal/OpenFOAM/harshal-.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/',100,'c',5,'SST')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/46_1Inlet10/postProcessing/',12,'r',5,'46.1$^o$ $\gamma Re_{\\theta}$')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/35Inlet10/postProcessing/',30,'b',5,'35$^o$ $\gamma Re_{\\theta}$')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/FinerDNSgrid/TransSteady/Tu4Suluksna/postProcessing/',10,'g',5,'Finer grid $\gamma Re_{\\theta}$','yes')

#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percent/postProcessing/',9.3,'b',3,'yPlus')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/30Inlet10/postProcessing/',5,'r',5,'30$^o$')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/40Inlet10/postProcessing/',30,'r',5,'40$^o$')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/43Inlet10/postProcessing/',30,'g',5,'43$^o$')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',38.8,'r',5,'LKE 38800')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',38.8,'m',9,'LKE 38800')
##plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',52.4,'brown',9,'LKE 52400')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',70,'grey',9,'LKE 70000')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinear/postProcessing/',0.022,'b--',9,'LKEURANS22')

##plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinear/postProcessing/',0.16,'m--',9,'LKEURANS160')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinear/postProcessing/',0.03,'c--',9,'LKEURANS30')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinear/postProcessing/',0.18,'r--',9,'LKEURANS180')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinear/postProcessing/',0.2,'c*',9,'LKEURANS200')

#plotWallShearStress('/home/harshal/Open6FOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/',100,'k',9,'SST')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.18,'y^',9,'Laminar18')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.34,'y',9,'Laminar34')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.66,'ro',9,'Laminar66')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.74,'bo',9,'Laminar74')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.9801,'y',9,'Laminar98')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.9001,'g',9,'Laminar90')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/URANS/Laminar/LaminarCaseFeb2107/postProcessing/',0.8201,'v',9,'Laminar82')
#
##plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.90/run/RAS/T106A/LPT_DNSGrid/TransitionURANS/LKE/Tu4perRANSLimLinearTry2/postProcessing/',0.022,'y*',9,'LKEURANS')
##plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/46_1Inlet10/poexstProcessing/',12,'r',9,'46.1$^o$ $\gamma Re_{\\theta}$')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/35Inlet10/postProcessing/',30,'b',9,'35$^o$ $\gamma Re_{\\theta}$')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/FinerDNSgrid/TransSteady/Tu4Suluksna/postProcessing/',10,'g',9,'Finer grid $\gamma Re_{\\theta}$','yes')


#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',8,'r',6,'SST scaled 8')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',74,'b',6,'SST')


#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentSpart/postProcessing/',2.188,'grey',9,'LKE 2.188 ')
#
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/',12.8,'g',9,'LKE 4.8 ')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/',6,'grey',9,'LKE 4.8 ')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',5.792,'r',9,'SST upwind ')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',74,'r',9,'SST (Turb)')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',7,'b',9, '$\gamma Re_{\\theta}$')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',14,'b',9, '$\gamma Re_{\\theta}$ (Trans)')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',0.04201,'m',9, 'URANS$\gamma Re_{\\theta} 2400$ (Trans)')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',0.04201,'m--',9, 'URANS$\gamma Re_{\\theta} 4201$ (Trans)')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',50,'y',9, ' Sul 50$\gamma Re_{\\theta}$')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/',50,'c',9, 'Langtry$\gamma Re_{\\theta}$')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/',6,'g',9, 'Laminar')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/',14,'brown',9, 'Laminar14')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_0001/postProcessing/',100,'grey',9, '0_0001')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_01/postProcessing/',40,'brown',9, '0_02',Re=60000) 
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/Nowallfunc_0_0001/postProcessing/',14,'r',9,'No wall func')
#100k cases
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/',8,'m--',9, 'SST 100k 8',Re=100000)
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/Tu4LangtryGradOff/postProcessing/',8,'b--',9,'100k $\gamma Re_{\\theta} 5$',Re=100000)
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/',100,'k',9,'LKE 100k',Re=100000)
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/URANS/laminar/postProcessing/',0.0501,'k*',9,'laminar 60k URANS av',Re=60000)
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/URANS/laminar/postProcessing/',0.05,'m*',9,'laminar 60k URANS',Re=60000)
#plotTauwData(9) 
#plotTauwDataRe100kTu4(9)



