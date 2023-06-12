# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:46:50 2016

@author: harshal

"""
import PostProcess
from numpy import genfromtxt
import matplotlib.pyplot as plt
plt.close("all")
Loading = 'Schemes'
wakeLoss = 'Schemes'
#timeConvergenceDirectory = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilonLimitedLinear/postProcessing/'
#timeConvergenceDirectory = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/'
#timeConvergence(timeConvergenceDirectory,'SST',time)
#timeConvergenceDirectory = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilon/postProcessing/'
#time = [30,60,90] 5
casesLoading(Loading)
casesWakeLoss(wakeLoss)

#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',52.4,'m',5,'LKE 52400')
#plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/',100,'c',5,'SST')
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
plotyPlus('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',74,'b',6,'SST')


#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentSpart/postProcessing/',2.188,'grey',9,'LKE 2.188 ')
#
plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/',4.8,'g',9,'LKE 4.8 ')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/LKE/Tu4percentLocal/postProcessing/',6,'grey',9,'LKE 4.8 ')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',5.792,'r',9,'SST upwind ')
plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/SST/Tu4percent/postProcessing/',74,'r',9,'SST (Turb)')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',7,'b',9, '$\gamma Re_{\\theta}$')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',14,'b',9, '$\gamma Re_{\\theta}$ (Trans)')
plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',0.024,'m',9, 'URANS$\gamma Re_{\\theta} 2400$ (Trans)')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',0.04201,'m--',9, 'URANS$\gamma Re_{\\theta} 4201$ (Trans)')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4percentSuluksna/postProcessing/',50,'y',9, ' Sul 50$\gamma Re_{\\theta}$')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/Tu4LangtryGradOff/postProcessing/',50,'c',9, 'Langtry$\gamma Re_{\\theta}$')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/',6,'g',9, 'Laminar')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/Steady/laminar/postProcessing/',14,'brown',9, 'Laminar14')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_0001/postProcessing/',100,'grey',9, '0_0001')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/0_01/postProcessing/',40,'brown',9, '0_02',Re=60000) 
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/TransSteady/GammaReTheta/VaryInletLengthScale/Nowallfunc_0_0001/postProcessing/',14,'r',9,'No wall func')
#100k cases
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/Steady/SST/Tu4percent/postProcessing/',8,'c--',9, 'SST 100k 8',Re=100000)
plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/Tu4LangtryGradOff/postProcessing/',8,'b--',9,'100k $\gamma Re_{\\theta} 5$',Re=100000)
plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPTDNSGridScaled/100k/TransSteady/LKE/Tu4percentLocal/postProcessing/',100,'r',9,'LKE 100k',Re=100000)
plotTauwData(9) 

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/',10.8,'m',9,'LKE 10800')

#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/',90,'r',89,'SST')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percent/postProcessing/',9.3,'m',89,'Transition')
#plotWallShearStress('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/VaryInletAngleLengthScale/35Inlet10/postProcessing/',30,'b',89,'TM')
#tcD = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/realizablekEpsilonLimitedLinear/postProcessing/'
#tcD = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/'
#timeConvergence(tcD,'SST',time)
#timeConvergence(tcD,'kEpsilon',time)

#f2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/sets/70/line2_U_wallShearStress.xy'
#velocity = np.loadtxt(f2)
#MagVelo = (velocity[:,1]**2 + velocity[:,2]**2)**0.5
#plt.figure(900)
#plt.plot(MagVelo/85.4,(velocity[:,0]-velocity[0,0])/0.5,'r',linewidth=2,label = 'RANS Lam')
#
#f3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percentSuluksna/postProcessing/sets/3.5/line2_U.xy'
#velocity = np.loadtxt(f3)
#MagVelo = (velocity[:,1]**2 + velocity[:,2]**2)**0.5
#plt.plot(MagVelo/85.4,(velocity[:,0]-velocity[0,0])/0.5,linewidth=2,label = 'RANS turb')
#a = [ -3.43947664e-10,   7.26856418e-03,   2.90701672e-02,
#         6.54099705e-02,   1.15577469e-01,   1.79119645e-01,
#         2.54570681e-01,   3.41374686e-01,   4.38197431e-01,
#         5.42793576e-01,   6.54017140e-01,   7.69200226e-01,
#         8.86406872e-01,   1.00232485e+00,   1.11534286e+00,
#         1.22343857e+00,   1.32482405e+00,   1.41733339e+00,
#         1.50007822e+00,   1.57145458e+00,   1.62967357e+00,
#         1.67461450e+00,   1.70674690e+00,   1.72724094e+00,
#         1.73880086e+00,   1.74380781e+00,   1.74433939e+00,
#         1.74197332e+00,   1.73778796e+00,   1.73246773e+00,
#         1.72642766e+00,   1.71991210e+00,   1.71306868e+00,
#         1.70597420e+00,   1.69866878e+00,   1.69117081e+00,
#         1.68349208e+00,   1.67564508e+00,   1.66763865e+00,
#         1.65947876e+00,   1.65116349e+00,   1.64267936e+00,
#         1.63402707e+00,   1.62523221e+00,   1.61634511e+00,
#         1.60739921e+00,   1.59839635e+00,   1.58932394e+00,
#         1.58016105e+00,   1.57090571e+00]
#         
#b = [0.00000000e+00,   2.08265572e-05,   8.32945301e-05,
#         1.87418740e-04,   3.33196060e-04,   5.20614798e-04,
#         7.49689772e-04,   1.02040617e-03,   1.33277567e-03,
#         1.68679828e-03,   2.08246544e-03,   2.51978570e-03,
#         2.99874739e-03,   3.51936218e-03,   4.08163322e-03,
#         4.68554567e-03,   5.33111435e-03,   6.01832133e-03,
#         6.74718454e-03,   7.51770087e-03,   8.32985861e-03,
#         9.18367259e-03,   1.00791280e-02,   1.10162365e-02,
#         1.19950012e-02,   1.30154043e-02,   1.40774667e-02,
#         1.51811790e-02,   1.63265360e-02,   1.75135312e-02,
#         1.87421826e-02,   2.00124903e-02,   2.13244511e-02,
#         2.26780534e-02,   2.40732970e-02,   2.55102086e-02,
#         2.69887615e-02,   2.85089559e-02,   3.00708035e-02,
#         3.16742924e-02,   3.33194524e-02,   3.50062506e-02,
#         3.67346903e-02,   3.85047863e-02,   4.03165353e-02,
#         4.21699289e-02,   4.40649725e-02,   4.60016754e-02,
#         4.79800018e-02,   4.99999993e-02]
#plt.plot(a,b,'k--',linewidth=2,label = 'DNS')
#plt.xlabel('${U}/{U_0}$')
#plt.ylabel('$y/y_0$')
#plt.legend(loc=0)
#
#
#f4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/LKEModel/Tu4percent/postProcessing/sets/70/line3_U_wallShearStress.xy'
#velocity = np.loadtxt(f4)
#MagVelo = (velocity[:,1]**2 + velocity[:,2]**2)**0.5
#plt.figure(902)
#plt.plot(MagVelo/85.4,(velocity[:,0]-velocity[0,0])/0.1,'r',linewidth=2,label = 'RANS Lam')
#
#f5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/TransitionSteady/Tu4percentSuluksna/postProcessing/sets/3.5/line3_U.xy'
#velocity = np.loadtxt(f5)
#MagVelo = (velocity[:,1]**2 + velocity[:,2]**2)**0.5
#plt.plot(MagVelo/85.4,(velocity[:,0]-velocity[0,0])/0.1,linewidth=2,label = 'RANS turb')
#
#a1 = [  3.92489267e-10,   8.09709915e-04,   3.23869156e-03,
#         7.28895185e-03,   1.34722798e-02,   2.21177403e-02,
#         3.43138760e-02,   5.04436721e-02,   7.14357590e-02,
#         9.87142954e-02,   1.32862216e-01,   1.75360327e-01,
#         2.27099777e-01,   2.89633910e-01,   3.63602513e-01,
#         4.49578705e-01,   5.47746483e-01,   6.57520322e-01,
#         7.77311195e-01,   9.04762761e-01,   1.03649800e+00,
#         1.16818129e+00,   1.29539875e+00,   1.41361171e+00,
#         1.51888876e+00,   1.60848185e+00,   1.68135403e+00,
#         1.73779366e+00,   1.77950432e+00,   1.80893019e+00,
#         1.82878485e+00,   1.84151346e+00,   1.84922710e+00,
#         1.85363516e+00,   1.85579787e+00,   1.85654694e+00,
#         1.85636180e+00,   1.85558097e+00,   1.85439233e+00,
#         1.85293755e+00,   1.85130237e+00,   1.84954659e+00,
#         1.84771582e+00,   1.84584305e+00,   1.84394716e+00,
#         1.84203362e+00,   1.84010448e+00,   1.83815763e+00,
#         1.83617435e+00,   1.83413367e+00]
#         
#b1 =  [0.00000000e+00,   2.08225494e-05,   8.32864895e-05,
#         1.87443526e-04,   3.33193955e-04,   5.20637479e-04,
#         7.49674398e-04,   1.02040441e-03,   1.33277952e-03,
#         1.68679231e-03,   2.08245391e-03,   2.51980490e-03,
#         2.99874928e-03,   3.51939046e-03,   4.08162133e-03,
#         4.68554529e-03,   5.33111435e-03,   6.01832481e-03,
#         6.74717665e-03,   7.51772529e-03,   8.32986362e-03,
#         9.18369876e-03,   1.00791236e-02,   1.10162415e-02,
#         1.19950045e-02,   1.30154089e-02,   1.40774547e-02,
#         1.51811973e-02,   1.63265333e-02,   1.75135624e-02,
#         1.87421811e-02,   2.00124930e-02,   2.13244500e-02,
#         2.26780483e-02,   2.40732918e-02,   2.55102283e-02,
#         2.69887583e-02,   2.85089813e-02,   3.00707940e-02,
#         3.16743036e-02,   3.33194545e-02,   3.50062468e-02,
#         3.67346842e-02,   3.85047630e-02,   4.03165386e-02,
#         4.21699556e-02,   4.40649660e-02,   4.60016676e-02,
#         4.79800143e-02,   5.00000024e-02 ]       
#
#plt.plot(a,b,'k--',linewidth=2,label = 'DNS')
#plt.xlabel('${U}/{U_0}$')
#plt.ylabel('$y/y_0$')
#plt.legend(loc=0)