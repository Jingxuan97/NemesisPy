from NemesisPy.Profile import *
from NemesisPy.Models import *
from NemesisPy.Data import *
import numpy as np
import matplotlib.pyplot as plt

###############################################################################################

def npvar(nvar,npro,varident,varparam):
    
    """
        FUNCTION NAME : npvar()
        
        DESCRIPTION :
        
            Function for getting the number of points associated with each variable in Nemesis
        
        INPUTS :
        
            nvar :: Number of variables to retrieve
            npro :: Number of altitude points in atmospheric profiles
            varident(nvar,3) :: Variable ID
            varparam(nvar,mparam) :: Extra parameters for describing the retrieved variable
        
        OPTIONAL INPUTS: none
        
        OUTPUTS :
        
            nxvar(nvar) :: Number of points in state vector associated with each variable
        
        CALLING SEQUENCE:
        
            nxvar = npvar(nvar,npro,varident,varparam)
        
        MODIFICATION HISTORY : Juan Alday (29/04/2019)
        
    """
    
    nxvar = np.zeros(nvar,dtype='int32')
    for i in range(nvar):
        #if nvar>1:
        imod = varident[i,2]
        ipar = varparam[i,0]
        
        if imod == -1:
            nxvar[i] = npro
        elif imod == 0:
            nxvar[i] = npro
        elif imod == 1:
            nxvar[i] = 2
        elif imod == 2:
            nxvar[i] = 1
        elif imod == 3:
            nxvar[i] = 1
        elif imod == 4:
            nxvar[i] = 3
        elif imod == 5:
            nxvar[i] = 1
        elif imod == 6:
            nxvar[i] = 2
        elif imod == 7:
            nxvar[i] = 2
        elif imod == 8:
            nxvar[i] = 3
        elif imod == 9:
            nxvar[i] = 3
        elif imod == 10:
            nxvar[i] = 4
        elif imod == 11:
            nxvar[i] = 2
        elif imod == 12:
            nxvar[i] = 3
        elif imod == 13:
            nxvar[i] = 3
        elif imod == 14:
            nxvar[i] = 3
        elif imod == 15:
            nxvar[i] = 3
        elif imod == 16:
            nxvar[i] = 4
        elif imod == 17:
            nxvar[i] = 2
        elif imod == 18:
            nxvar[i] = 2
        elif imod == 19:
            nxvar[i] = 4
        elif imod == 20:
            nxvar[i] = 2
        elif imod == 21:
            nxvar[i] = 2
        elif imod == 22:
            nxvar[i] = 5
        elif imod == 23:
            nxvar[i] = 4
        elif imod == 24:
            nxvar[i] = 3
        elif imod == 25:
            nxvar[i] = int(ipar)
        elif imod == 26:
            nxvar[i] = 4
        elif imod == 27:
            nxvar[i] = 3
        elif imod == 28:
            nxvar[i] = 1
        elif imod == 228:
            nxvar[i] = 7
        elif imod == 229:
            nxvar[i] = 7
        elif imod == 230:
            nxvar[i] = 2*int(ipar)
        elif imod == 444:
            nxvar[i] = 1 + 1 + int(ipar)
        elif imod == 666:
            nxvar[i] = 1
        elif imod == 999:
            nxvar[i] = 1
        else:
            sys.exit('error :: varID not included in npvar_nemesis()')

    return nxvar

###############################################################################################

def modelm1(atm,ipar,xprof,MakePlot=False):
    
    """
        FUNCTION NAME : model0()
        
        DESCRIPTION :
        
            Function defining the model parameterisation -1 in NEMESIS.
            In this model, the aerosol profiles is modelled as a continuous profile in units
            of particles per cm3. Note that typical units of aerosol profiles in NEMESIS
            are in particles per gram of atmosphere
        
        INPUTS :
        
            atm :: Python class defining the atmosphere

            ipar :: Atmospheric parameter to be changed
                    (0 to NVMR-1) :: Gas VMR
                    (NVMR) :: Temperature
                    (NVMR+1 to NVMR+NDUST-1) :: Aerosol density
                    (NVMR+NDUST) :: Para-H2
                    (NVMR+NDUST+1) :: Fractional cloud coverage

            xprof(npro) :: Atmospheric aerosol profile in particles/cm3
        
        OPTIONAL INPUTS:

            MakePlot :: If True, a summary plot is generated
        
        OUTPUTS :
        
            atm :: Updated atmospheric class
            xmap(npro,ngas+2+ncont,npro) :: Matrix of relating funtional derivatives to 
                                             elements in state vector
        
        CALLING SEQUENCE:
        
            atm,xmap = modelm1(atm,ipar,xprof)
        
        MODIFICATION HISTORY : Juan Alday (29/03/2021)
        
    """

    npro = len(xprof)
    if npro!=atm.NP:
        sys.exit('error in model 0 :: Number of levels in atmosphere does not match and profile')

    npar = atm.NVMR+2+atm.NDUST
    xmap = np.zeros([npro,npar,npro])

    if ipar<atm.NVMR:  #Gas VMR
        sys.exit('error :: Model -1 is just compatible with aerosol populations')
    elif ipar==atm.NVMR: #Temperature
        sys.exit('error :: Model -1 is just compatible with aerosol populations')
    elif ipar>atm.NVMR:
        jtmp = ipar - (atm.NVMR+1)
        x1 = np.exp(xprof)
        if jtmp<atm.NDUST:
            rho = atm.calc_rho(molwt)  #kg/m3
            rho = rho / 1.0e3 #g/cm3
            atm.DUST[:,jtmp] = x1 / rho
        elif jtmp==atm.NDUST:
            sys.exit('error :: Model -1 is just compatible with aerosol populations')
        elif jtmp==atm.NDUST+1:
            sys.exit('error :: Model -1 is just compatible with aerosol populations')
    
    for j in range(npro):
        xmap[0:npro,ipar,j] = x1[:] / rho[:]
        

    if MakePlot==True:
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

        for i in range(atm.NDUST):
            ax1.semilogx(atm.DUST[:,i]*rho,atm.H)
            ax2.semilogx(atm.DUST[:,i],atm.H)

        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax1.set_xlabel('Aerosol density (particles per cm$^{-3}$)')
        ax1.set_ylabel('Altitude (km)')
        ax2.set_xlabel('Aerosol density (particles per gram of atm)')
        ax2.set_ylabel('Altitude (km)')
        plt.tight_layout()
        plt.show()

    return atm,xmap


###############################################################################################

def model0(atm,ipar,xprof,MakePlot=False):
    
    """
        FUNCTION NAME : model0()
        
        DESCRIPTION :
        
            Function defining the model parameterisation 0 in NEMESIS.
            In this model, the atmospheric parameters are modelled as continuous profiles
            in which each element of the state vector corresponds to the atmospheric profile 
            at each altitude level
        
        INPUTS :
        
            atm :: Python class defining the atmosphere

            ipar :: Atmospheric parameter to be changed
                    (0 to NVMR-1) :: Gas VMR
                    (NVMR) :: Temperature
                    (NVMR+1 to NVMR+NDUST-1) :: Aerosol density
                    (NVMR+NDUST) :: Para-H2
                    (NVMR+NDUST+1) :: Fractional cloud coverage

            xprof(npro) :: Atmospheric profile
        
        OPTIONAL INPUTS:

            MakePlot :: If True, a summary plot is generated
        
        OUTPUTS :
        
            atm :: Updated atmospheric class
            xmap(npro,ngas+2+ncont,npro) :: Matrix of relating funtional derivatives to 
                                             elements in state vector
        
        CALLING SEQUENCE:
        
            atm,xmap = model0(atm,ipar,xprof)
        
        MODIFICATION HISTORY : Juan Alday (29/03/2021)
        
    """

    npro = len(xprof)
    if npro!=atm.NP:
        sys.exit('error in model 0 :: Number of levels in atmosphere does not match and profile')

    npar = atm.NVMR+2+atm.NDUST
    xmap = np.zeros([npro,npar,npro])

    if ipar<atm.NVMR:  #Gas VMR
        jvmr = ipar
        x1 = np.exp(xprof)
        atm.VMR[:,jvmr] = x1
    elif ipar==atm.NVMR: #Temperature
        x1 = xprof
        atm.T[:] = x1
    elif ipar>atm.NVMR:
        jtmp = ipar - (atm.NVMR+1)
        x1 = np.exp(xprof)
        if jtmp<atm.NDUST:
            atm.DUST[:,jtmp] = x1
        elif jtmp==atm.NDUST:
            atm.PARAH2 = x1
        elif jtmp==atm.NDUST+1:
            atm.FRAC = x1
    
    for j in range(npro):
        xmap[0:npro,ipar,j] = x1[:]
        

    if MakePlot==True:
        fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,5))

        ax1.semilogx(atm.P/101325.,atm.H)
        ax2.plot(atm.T,atm.H)
        for i in range(atm.NVMR):
            ax3.semilogx(atm.VMR[:,i],atm.H)

        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax1.set_xlabel('Pressure (atm)')
        ax1.set_ylabel('Altitude (km)')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Altitude (km)')
        ax3.set_xlabel('Volume mixing ratio')
        ax3.set_ylabel('Altitude (km)')
        plt.tight_layout()
        plt.show()

    return atm,xmap


###############################################################################################

def model2(atm,ipar,scf,MakePlot=False):
    
    """
        FUNCTION NAME : model2()
        
        DESCRIPTION :
        
            Function defining the model parameterisation 2 in NEMESIS.
            In this model, the atmospheric parameters are scaled using a single factor with 
            respect to the vertical profiles in the reference atmosphere
        
        INPUTS :
        
            atm :: Python class defining the atmosphere

            ipar :: Atmospheric parameter to be changed
                    (0 to NVMR-1) :: Gas VMR
                    (NVMR) :: Temperature
                    (NVMR+1 to NVMR+NDUST-1) :: Aerosol density
                    (NVMR+NDUST) :: Para-H2
                    (NVMR+NDUST+1) :: Fractional cloud coverage

            scf :: Scaling factor
        
        OPTIONAL INPUTS:

            MakePlot :: If True, a summary plot is generated
        
        OUTPUTS :
        
            atm :: Updated atmospheric class
            xmap(1,ngas+2+ncont,npro) :: Matrix of relating funtional derivatives to 
                                             elements in state vector
        
        CALLING SEQUENCE:
        
            atm,xmap = model2(atm,ipar,scf)
        
        MODIFICATION HISTORY : Juan Alday (29/03/2021)
        
    """

    npar = atm.NVMR+2+atm.NDUST
    xmap = np.zeros([1,npar,atm.NP])

    x1 = np.zeros(atm.NP)
    if ipar<atm.NVMR:  #Gas VMR
        jvmr = ipar
        x1[:] = atm.VMR[:,jvmr] * scf
        atm.VMR[:,jvmr] =  x1
    elif ipar==atm.NVMR: #Temperature
        x1[:] = atm.T[:] * scf
        atm.T[:] = x1 
    elif ipar>atm.NVMR:
        jtmp = ipar - (atm.NVMR+1)
        if jtmp<atm.NDUST:
            x1[:] = atm.DUST[:,jtmp] * scf
            atm.DUST[:,jtmp] = x1
        elif jtmp==atm.NDUST:
            x1[:] = atm.PARAH2 * scf
            atm.PARAH2 = x1
        elif jtmp==atm.NDUST+1:
            x1[:] = atm.FRAC * scf
            atm.FRAC = x1

    xmap[0,ipar,:] = x1[:]
    
    if MakePlot==True:
        fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,5))

        ax1.semilogx(atm.P/101325.,atm.H)
        ax2.plot(atm.T,atm.H)
        for i in range(atm.NVMR):
            ax3.semilogx(atm.VMR[:,i],atm.H)

        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax1.set_xlabel('Pressure (atm)')
        ax1.set_ylabel('Altitude (km)')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Altitude (km)')
        ax3.set_xlabel('Volume mixing ratio')
        ax3.set_ylabel('Altitude (km)')
        plt.tight_layout()
        plt.show()

    return atm,xmap


###############################################################################################

def model3(atm,ipar,scf,MakePlot=False):
    
    """
        FUNCTION NAME : model2()
        
        DESCRIPTION :
        
            Function defining the model parameterisation 2 in NEMESIS.
            In this model, the atmospheric parameters are scaled using a single factor 
            in logscale with respect to the vertical profiles in the reference atmosphere
        
        INPUTS :
        
            atm :: Python class defining the atmosphere

            ipar :: Atmospheric parameter to be changed
                    (0 to NVMR-1) :: Gas VMR
                    (NVMR) :: Temperature
                    (NVMR+1 to NVMR+NDUST-1) :: Aerosol density
                    (NVMR+NDUST) :: Para-H2
                    (NVMR+NDUST+1) :: Fractional cloud coverage

            scf :: Log scaling factor
        
        OPTIONAL INPUTS:

            MakePlot :: If True, a summary plot is generated
        
        OUTPUTS :
        
            atm :: Updated atmospheric class
            xmap(1,ngas+2+ncont,npro) :: Matrix of relating funtional derivatives to 
                                             elements in state vector
        
        CALLING SEQUENCE:
        
            atm,xmap = model2(atm,ipar,scf)
        
        MODIFICATION HISTORY : Juan Alday (29/03/2021)
        
    """

    npar = atm.NVMR+2+atm.NDUST
    xmap = np.zeros([1,npar,atm.NP])

    x1 = np.zeros(atm.NP)
    if ipar<atm.NVMR:  #Gas VMR
        jvmr = ipar
        x1[:] = atm.VMR[:,jvmr] * np.exp(scf)
        atm.VMR[:,jvmr] =  x1 
    elif ipar==atm.NVMR: #Temperature
        x1[:] = atm.T[:] * np.exp(scf)
        atm.T[:] = x1 
    elif ipar>atm.NVMR:
        jtmp = ipar - (atm.NVMR+1)
        if jtmp<atm.NDUST:
            x1[:] = atm.DUST[:,jtmp] * np.exp(scf)
            atm.DUST[:,jtmp] = x1
        elif jtmp==atm.NDUST:
            x1[:] = atm.PARAH2 * np.exp(scf)
            atm.PARAH2 = x1
        elif jtmp==atm.NDUST+1:
            x1[:] = atm.FRAC * np.exp(scf)
            atm.FRAC = x1

    xmap[0,ipar,:] = x1[:]
    
    if MakePlot==True:
        fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,5))

        ax1.semilogx(atm.P/101325.,atm.H)
        ax2.plot(atm.T,atm.H)
        for i in range(atm.NVMR):
            ax3.semilogx(atm.VMR[:,i],atm.H)

        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax1.set_xlabel('Pressure (atm)')
        ax1.set_ylabel('Altitude (km)')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Altitude (km)')
        ax3.set_xlabel('Volume mixing ratio')
        ax3.set_ylabel('Altitude (km)')
        plt.tight_layout()
        plt.show()

    return atm,xmap