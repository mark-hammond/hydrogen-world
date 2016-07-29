#Modified for Exoclimes2012.  Uses dry enthalpy conserving convective
#adjustment and implements an energy conserving surface flux scheme
#(Note enthalpy is the right thing to conserve, not dry static energy)
#
#Modified further to make it easier to swap in homebrew or graygas models.
# (8/28/2012, for Beijing lectures)
#
#
#Note: The stellar flux past the shortwave cutoff had to be
#added in to the calculation of incoming flux. This needs
#to be put in the main code-tree.

#ToDo:
#     *Make it possible to use ccm radiation.
#         ->Need radcomp in ccmradFunctions. Return heating in W/kg
#     *make it possible to use the fancy version of miniclimt
#     *Make a script that computes the pure radiative equilibria
#      in Figure {fig:RealGasRadEq}, and also the logarithmic slopes.
#     *Handle the contribution of water vapor to surface pressure
#      in setting up the pressure grid.  (Tricky, because
#      temperature is changing!)  That's important when
#      surface temperatures are much over 300K.  The changing
#      surface pressure makes the time stepping rather tricky.

#Data on section of text which this script is associated with
Chapter = '5'
Section = '**'
Figure = '**'
#

#import climt_lite as climt #ccm radiation model (Replace with ccmradFunctions?)
import GreyHeat as Grey
from ClimateUtilities import *
import math,phys
import planets
import matplotlib.pyplot as plt


#Dry adjustment routine.
#**ToDo: Modify so it handles non-uniform pressure levels
#        correctly, and conserves integrated enthalpy cp T
#**ToDo: Update this so that it conserves dry static energy.
#
#Iterative routine for dry convective adjustment
def dryAdj(T,p):
    #Downward pass
    for i in range(len(T)-1):
        if doVarRcp:
            Rcp = 8.31 / shomateRef(T[i])
        else:
            Rcp = 2./7.
        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        pfact = (p1/p2)**Rcp
        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) #Equal layer masses
                              #Not quite compatible with how
                              #heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            T[i] = T1
            T[i+1] = T2
    #Upward pass
    for i in range(len(T)-2,-1,-1):
        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        pfact = (p1/p2)**Rcp
        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) #Equal layer masses
                              #Not quite compatible with how
                              #heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            T[i] = T1
            T[i+1] = T2
            
#Shomate equation to give Rcp in different T ranges
def shomate(temp):
    if temp < 1000:
        A, B, C, D, E = 33.066178, -11.363417, 11.432816, -2.772874, -0.158558
    elif temp < 2500:
        A, B, C, D, E = 18.563083, 12.257357, -2.859786, 0.268238, 1.977990	 
    else:
        A, B, C, D, E = 43.413560, -4.293079, 1.272428, -0.096876, -20.533862
    t = temp / 1000.
    Rcp = A + B*t + C*(t**2) + (D*t**3) + E/(t**2)
    return Rcp
    

    
#--------------------Set radmodel options-------------------
#---Instantiate the radiation model---
#r = climt.radiation()
#n = r.nlev
n = 20

#Set global constants
ps = 1000.
rh = 1.e-30#Relative humidity
rhbdd = 1.e-30
dt = 24.*3600. #time step in seconds

#---Set up pressure array (a global)----
ptop = 50. #Top pressure in mb (Changed from 1mb in original)
pstart = .995*ps
rat = (ptop/pstart)**(1./n)
logLevels = [pstart*rat**i for i in range(n)]
logLevels.reverse()
levels = [ptop + i*(pstart-ptop)/(n-1) for i in range(n)]
#p=Numeric.array(levels,Numeric.Float)
p = numpy.array(logLevels)

#---Temperature and moisture arrays (initialized)
#  **CHANGE: Eliminated typecode for numpy compatibility
T = numpy.zeros(n) + 230.
q = numpy.zeros(n) + 1.e-30
small = numpy.zeros(n) + 1.e-30


###Set the ozone profile
###r(p,ps,T,300.,q)
###o3 = r.o3
##o3max = 2.e-5 #Typical value is 2.e-5
### **CHANGE: Eliminated typecode
##o3= Numeric.zeros(r.nlev) + 1.e-30
##for i in range(n):
##    o3[i] = max(o3max*math.exp(-(p[i]-20.)**2/50.**2),1.e-20)

#CCM radiative flux and heating. 
def radcompCCM(p,T,Tg,q):
    #**Changed call to r
    r(p=p,ps=ps,T=T,Ts=Tg,q=q,co2 = co2,o3=small)
    #**Changed indexing of flux and hr
    fluxLW,heatLW = -r.lwflx[:,0,0],r.lwhr[:,0,0]
    fluxStellar,heatStellar = -r.swflx[:,0,0],r.swhr[:,0,0]
    flux = fluxLW+fluxStellar
    if doStellarAbs:    
        heat = heatLW+heatStellar
    else:
        heat = heatLW
    return flux,heat

#Define function to do time integration for n steps
def steps(Tg,T,nSteps,dtime):
    for i in range(nSteps):
        #Tg = T[-1] #let it float
        #Do smoothing
##        if i%20 == 0:
##            for j in range(1,len(T)-1):
##                T[j] = .25*T[j-1] + .5*T[j] + .25*T[j+1]
        flux,heat = Grey.radcompLW(p,T,Tg,q) #radcompCCM(p,T,Tg,q)
        dT = heat*dtime
        #Limit the temperature change per step
        dT = numpy.where(dT>5.,5.,dT)
        dT = numpy.where(dT<-5.,-5.,dT)
        #Midpoint method time stepping
        #changed call to r.  Also modified to hold Tg fixed
        fluxLW,heatLW = Grey.radcompLW(p,T+.5*dT,Tg,q)
        fluxStellar,heatStellar = Grey.radcompStellar(p,T+.5*dT,Tg,q)
        if doStellarAbs:
            heat = heatLW + heatStellar
            flux = fluxLW + fluxStellar
        else:
            flux = fluxLW
            heat = heatLW
        dT = heat*dtime
        #Limit the temperature change per step
        dT = Numeric.where(dT>5.,5.,dT)
        dT = Numeric.where(dT<-5.,-5.,dT)
        T += dT
        #
        dTmax = max(abs(dT)) #To keep track of convergence
        
        #**Replace following with surface balance
##        
##        #Using flux[0] in the following is not a mistake.
##        #This is just a trick to relax towards satisfying
##        #the TOA balance. The formulation assumes that
##        #turbulent fluxes keep the low level air temperature
##        #near the ground temperature.
##
##        #Estimate the surface temperature change needed to
##        #bring the TOA energy budget into balance
##        Trad = (fluxLW[0]/phys.sigma)**.25
##        dTg = -flux[0]/(4.*phys.sigma*Trad**3)
##        #Relax partway towards that value
##        Tg = Tg + .1*dTg
        #
        #Replace following with pairwise convective adjustment
        #Tad = Tg*(p/p[-1])**Rcp
        #T = Numeric.where(T<Tad,Tad,T)
        #T[-1] = Tg
        #
        kturb = .1
        T[-1] += -dtime*kturb*(T[-1] - Tg)
        #Dry adjustment step
        for iadj in range(10):
            dryAdj(T,p)
        Tad = T[-1]*(p/p[-1])**Rcp
        #** Temporary kludge to keep stratosphere from getting too cold
        T = numpy.where(T<50.,50.,T)  #**KLUDGE
        #
        #Dummies for separate LW and stellar. **FIX THIS**
        #fluxStellar = fluxLW = heatStellar = heatLW = numpy.zeros(n)
    return Tg,Tad,T,flux,fluxStellar,fluxLW,heat,heatStellar,heatLW
    
#--------------Initializations-------------------------------------------
    
#----------------Set initial time step--------------------------------
#** Has something changed with climt? this seems to need smaller time step
#   I think maybe heating rate in climt is now K/day, not K/sec
dtime = 1# 1. # (for CO2 case; gray gas evolves faster)
            #Timestep in days; 5 days is the usual for Earthlike case
            #For the radiative convective case, you can get away with
            #using 50 for the first few hundred time steps, then
            #reducing to 5 or less later as the solution converges

#dtime = dtime*24.*3600. (Left over from old version of climt
#----------------------------------------------------------------------

#--------------------Set radmodel options-------------------
#---Instantiate the radiation model---
#r = climt.radiation()
#n = r.nlev
n = 20

#Set global constants
ps = 1000.
rh = 1.e-30#Relative humidity
rhbdd = 1.e-30
dt = 24.*3600. #time step in seconds
#---Set up pressure array (a global)----
ptop = 50. #Top pressure in mb (Changed from 1mb in original)
pstart = .995*ps
rat = (ptop/pstart)**(1./n)
logLevels = [pstart*rat**i for i in range(n)]
logLevels.reverse()
levels = [ptop + i*(pstart-ptop)/(n-1) for i in range(n)]
#p=Numeric.array(levels,Numeric.Float)
p = numpy.array(logLevels)


#--------------Other parameters-------------------------------------------
doStellarAbs = False
doVarRcp = False
#Set the gravity and thermodynamic constants
Rcp = 2./7.
#Set composition constants (globals)
co2 = 300.

#
#Initialize the temperature
#c = readTable(outputDir+'TTestTrop80.txt')
#p = c['p']
#T = c['T']
#Tg = T[-1]
Tg = 280.
#T = Tg*(p/p[-1])**Rcp
T = Tg*numpy.ones(len(p))
#T[-1] = Tg
#

#Set composition parameters for the radiation code you are using

#CCMrad:
# **CHANGE (climt): This is specified in the call to r
#r.params.value['co2'] = 300.

#Grey Gas:
Grey.tauInf = 5.
