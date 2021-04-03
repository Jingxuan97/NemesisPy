from NemesisPy.Profile import *
from NemesisPy.Models.Models import *
from NemesisPy.Data import *
import numpy as np
import matplotlib.pyplot as plt

###############################################################################################

def subprofretg(runname,atm,ispace,iscat,xlat,xlon,nvar,varident,varparam,nx,xn,jpre,flagh2p):

    """
    FUNCTION NAME : subprogretg()

    DESCRIPTION : Updates the atmosphere based on the variables and parameterisation in the
                  state vector

    INPUTS :
    
        runname :: Name of the Nemesis run
        atm :: Python class defining the atmosphere
        ispace :: (0) Wavenumber in cm-1 (1) Wavelength in um
        iscat :: Type of scattering calculation
        xlat :: Latitude of spectrum to be simulated
        xlon :: Longitude of spectrum to be simulated
        nvar :: Number of variables 
        varident(nvar,3) :: Variable or parameterisation ID 
        varparam(nvar,nparam) :: Additional parameters required by the parameterisation
        nx :: Number of elements in the state vector
        xn(nx) :: State vector 
        jpre :: Level of tanget pressure for limb pressure retrieval
        flagh2p :: Flag indicating whether para-H2 profile is variable

    OPTIONAL INPUTS: none
            
    OUTPUTS : 

        xmap(maxv,ngas+2+ncont,npro) :: Matrix relating functional derivatives calculated 
                                         by CIRSRADG to the elements of the state vector.
                                         Elements of XMAP are the rate of change of 
                                         the profile vectors (i.e. temperature, vmr prf
                                         files) with respect to the change in the state
                                         vector elements. So if X1(J) is the modified 
                                         temperature,vmr,clouds at level J to be 
                                         written out to runname.prf or aerosol.prf then
                                        XMAP(K,L,J) is d(X1(J))/d(XN(K)) and where
                                        L is the identifier (1 to NGAS+1+2*NCONT)

    CALLING SEQUENCE:

        xmap = subprofretg(runname,atm,ispace,iscat,xlat,xlon,nvar,varident,varparam,nx,xn,jpre,flagh2p)
 
    MODIFICATION HISTORY : Juan Alday (15/03/2021)

    """

    #Modify profile via hydrostatic equation to make sure the atm is in hydrostatic equilibrium
    if jpre!=-1:
        #Then we modify the altitude levels and keep the pressures fixed
        print('Finish this part')

    else:
        #Then we modifify the pressure levels and keep the altitudes fixed
        print('Finish this part')


    #Calculating molecular weight
    if atm.AMFORM==0:
        #We should pass the molecular weight from the .ref or .prf file here
        print('Finish this part')
    if atm.AMFORM==1:
        atm.adjust_VMR()
        molwt = atm.calc_molwt()
    elif atm.AMFORM==2:
        molwt = atm.calc_molwt()

    #Calculate atmospheric density
    rho = atm.calc_rho(molwt) #kg/m3

    #Look for the number of points in state vector associated with each variable
    nxvar = npvar(nvar,atm.NP,varident,varparam)
    
    #Initialising xmap
    xmap = np.zeros([nx,atm.NVMR+2+atm.NDUST,atm.NP])

    #Going through the different variables an updating the atmosphere accordingly
    ix = 0
    for ivar in range(nvar):

        if varident[ivar,0]<=100:
            
            #Reading the atmospheric profile which is going to be changed by the current variable
            xref = np.zeros([atm.NP])

            if varident[ivar,0]==0:     #Temperature is to be retrieved
                xref[:] = atm.T
                ipar = atm.NVMR
            elif varident[ivar,0]>0:    #Gas VMR is to be retrieved
                jvmr = np.where( (np.array(atm.ID)==varident[ivar,0]) & (np.array(atm.ISO)==varident[ivar,1]) )
                jvmr = int(jvmr[0])
                xref[:] = atm.VMR[:,jvmr]
                ipar = jvmr
            elif varident[ivar,0]<0:  
                jcont = int(-varident[ivar,0])
                if jcont>atm.NDUST+1:
                    sys.exit('error :: Variable outside limits',varident[ivar,0],varident[ivar,1],varident[ivar,2])
                elif jcont==atm.NDUST:   #Para-H2
                    if flagh2p==True:
                        xref[:] = atm.PARAH2
                    else:
                        sys.exit('error :: Para-H2 is declared as variable but atmosphere is not from Giant Planet')
                elif jcont==atm.NDUST+1: #Fractional cloud cover
                    xref[:] = atm.FRAC
                else:
                    xref[:] = atm.DUST[:,jcont]

                ipar = atm.NVMR + jcont


        x1 = np.zeros(atm.NP)        

        if varident[ivar,2]==-1:
#       Model -1. Continuous aerosol profile in particles cm-3
#       ***************************************************************

            xprof = np.zeros(nxvar[ivar])
            xprof[:] = xn[ix:ix+nxvar[ix]]
            atm,xmap1 = modelm1(atm,ipar,xprof,MakePlot=True)
            xmap[ix:ix+nxvar[ivar],:,0:atm.NP] = xmap1[:,:,:]

            ix = ix + nxvar[ivar]

        if varident[ivar,2]==0:
#       Model 0. Continuous profile
#       ***************************************************************

            xprof = np.zeros(nxvar[ivar])
            xprof[:] = xn[ix:ix+nxvar[ix]]
            atm,xmap1 = model0(atm,ipar,xprof,MakePlot=True)
            xmap[ix:ix+nxvar[ivar],:,0:atm.NP] = xmap1[:,:,:]

            ix = ix + nxvar[ivar]

        elif varident[ivar,2]==2:
#       Model 2. Scaling factor
#       ***************************************************************

            atm,xmap1 = model2(atm,ipar,xn[ix],MakePlot=True)
            xmap[ix:ix+nxvar[ivar],:,0:atm.NP] = xmap1[:,:,:]

            ix = ix + nxvar[ivar]

        elif varident[ivar,2]==3:
#       Model 2. Log scaling factor
#       ***************************************************************

            atm,xmap1 = model3(atm,ipar,xn[ix],MakePlot=True)
            xmap[ix:ix+nxvar[ivar],:,0:atm.NP] = xmap1[:,:,:]

            ix = ix + nxvar[ivar]

        else:
            sys.exit('error :: Model parameterisation has not yet been included')
