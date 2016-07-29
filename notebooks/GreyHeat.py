#===================================================================
#Computes fluxes and heating rates for the grey gas model.
#The fluxes are computed as a function of p/ps, given net optical
#thickness of the atmosphere tauinf .
#
#Modified to return dimensional heat and flux for use with
#the time-stepping radiative-convective model
#===================================================================

#Data on section of text which this script is associated with
#MODIFIED
Chapter = '4.**'
Figure = 'fig:AllTropNetIRFluxGrey'
#
#This is also the solution script for
Problem = '{Workbook:RadBalance2:PressBroadenedHeating}'
#This script can also be modified to use for the problem
# '{Workbook:RadBalance2:StratTropOLRGrey}'


import math,phys
from ClimateUtilities import *


#Grey gas transmission function.
#tauinf is a global
def Trans(tau1,tau2):
    return math.exp(-abs(tau1-tau2))

#Integrand for upward or downward flux. Note that
#the Schwartzschild integral is written here as an integral
#over p/ps, and correspondingly the gradient of T is written as
#dT/d(p/ps). The solution is written in the form of
#Eq. (4.13) (in First Edition).
#

def integrand(ppsp,params):
    #Without pressure broadening
    if PressureBroadening:
        tau1 = tauInf*(1.-ppsp**2)
        tau2 = tauInf*(1.-params.pps**2)
    else:
        tau1 = tauInf*(1.-ppsp)
        tau2 = tauInf*(1. - params.pps)
    Tfun = params.Tfun
    dTdp = (Tfun(ppsp+.01)-Tfun(ppsp-.01))/.02
    return Trans(tau1,tau2)*4.*phys.sigma*Tfun(ppsp)**3*dTdp

def Iplus(pps,Tfun,Tg):
    params = Dummy()
    params.pps = pps
    params.Tfun = Tfun
    Ts = Tfun(1.)
    limit = min(1.,pps+10./tauInf)
    quad = romberg(integrand,10)
    if PressureBroadening:
        tau = tauInf*(1.-pps**2)
    else:
        tau = tauInf*(1.-pps)
    BddTerm = (phys.sigma*Tg**4 - phys.sigma*Ts**4)*Trans(0.,tau)
    return quad([pps,limit],params,.1)+ phys.sigma*Tfun(pps)**4 +BddTerm

def Iminus(pps,Tfun,Tg):
    params = Dummy()
    params.pps = pps
    params.Tfun = Tfun
    limit = max(0.,pps-10./tauInf)
    quad = romberg(integrand,10)
    Tstrat = Tfun(0.)
    if PressureBroadening:
        tau = tauInf*(1.-pps**2)
    else:
        tau = tauInf*(1.-pps)
    return quad([pps,0.],params,.1)+ phys.sigma*Tfun(pps)**4 - phys.sigma*Tstrat**4*Trans(tau,tauInf)

#Return dimensional flux and heating
#pList is dimensional pressure.
#Moisture q not used for gray gas
#

#Pressure in mb, for consistency with ccmrad
def radcompLW(pList,TList,Tg,q):
    ppsL = pList/pList[-1] #p/psurf
    Tfun = interp(ppsL,TList)
    Ip = numpy.array([Iplus(pps,Tfun,Tg) for pps in ppsL])
    Im = numpy.array([Iminus(pps,Tfun,Tg) for pps in ppsL])
    flux = Ip-Im
    heat = -2.*phys.sigma*TList**4 + (Ip+Im)
    #
    #Re-dimensionalize heating to K/day
    if PressureBroadening :
        kappa = 2.*tauInf*g/(100.*pList[-1])
    else:
        kappa = tauInf*g/(100.*pList[-1])
    heat = heat*(kappa/cp)*24.*3600. #K/day
    return flux,heat

def radcompStellar(p,T,Tg,q,):
    n = len(p)
    ps = p[-1]
    Delta = 100000.
    wave = 2e4
    kappa = 0.05
    #
    #Mass ratio for computing molar concentration
    #  (used to compute proportion of self-collisions)
    massrat = GHG.MolecularWeight/BackgroundGas.MolecularWeight

    #Now compute the transmission function between each
    #level and the top.
    OrbitFact = (Lstellar/4.)/(phys.sigma*Tstellar**4)
    swFlux = OrbitFact*Delta*Planck(wave,Tstellar)
    flux = Numeric.zeros(n,Numeric.Float)
    flux[0] = -swFlux #Downward flux is negative
    #
    path = 0.
    pathq= 0.
    for j in range(1,n):
        moleCon = q[j]/(q[j] + (1-q[j])*massrat)
        #Temperature scaling------------
        #This is simplified to deal with the all-continuum
        #case
        wTS = (Tref/T[j])
        #---------------------------------
        #**Fix this; move to midpoint
        #path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
        path = path + wTS*(p[j]/pref)*abs(p[j-1]-p[j])/(g*cosThetaBar)
        #Compute the self-broadened pressure path 
        #pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
        Trans1 = TransStellar(path,kappa)
        flux[j] = -swFlux*Trans1
##########################################################################################
    #Flux computed. Now compute the heating
    heat = Numeric.zeros(n,Numeric.Float)
    #
    #Compute heat by taking the gradient of the flux.
    #The second order centered difference below allows
    #the computation to be done accurately even if
    #the pressure level spacing is non-uniform.
    #
    #Returns heating rate as W/kg.
    #Divide by specific heat to get K/s
    #
    delPlus = p[2:] - p[1:-1] #Should wind up as a [1:-1] array
    delMinus = p[1:-1]-p[:-2] #likewise
    A = (delMinus/delPlus)/(delPlus+delMinus)
    B = 1./delMinus - 1./delPlus
    C= (delPlus/delMinus)/(delPlus+delMinus)
    heat[1:-1] = A*flux[2:] + B*flux[1:-1] - C*flux[:-2]
    heat = g*heat #Convert to W/kg
    heat[0]=heat[1]
    heat[-1] = heat[-2]
    return flux,heat
    
def TransStellar(path,kappa):
    #kappa is the reference value at (pref,Tref)
    return math.exp(-path*kappa)
    
#Planck functions in wavenumber space (cm**-1)
def Planck(wavenum,T):
    return 100.*math.pi*phys.c*phys.B(100.*wavenum*phys.c,T)

#These are all globals
tauInf = 1.
PressureBroadening = True
GHG =phys.H2 #The greenhouse gas
BackgroundGas = phys.H2#The transparent background gas
pref = 1.e4 #Reference pressure at which band data is given
Tref = 260. #Reference temperature at which band data is given

#Set the photosphere temperature of the star
#(for shortwave radiation calculation
Tstellar = 3500.
#Set the stellar constant at the orbit
Lstellar = 4.*700.

#Gravity and cp (for diminensionalization)
g = 10.
cp = 1000.

#Mean slant path
cosThetaBar = .5









