# -*- coding: utf-8 -*-
import numpy as np
import Reference
from Radtrans import IMOD3



def GET_K(FILENAMES, NG, VWAVE, PRESS, TEMP, ID=None, ISO=None, IVAWE=None):
    """
    Inputs
    ------
    FILENAMES:
    PRESS:
    TEMP:
    IWAVE:
    VWAVE:

    Outputs
    -------
    KOUT(NWAVE, NLAYER, NGAS, NG):
        K-coefficients for each layer and each gas
    """
    NGAS = len(FILENAMES)
    NPOINTS = len(PRESS)
    NLAYER = NPOINTS
    NWAVE = len(VWAVE)
    assert len(FILENAMES) == NGAS
    KOUT = np.zeros((NWAVE, NLAYER, NGAS, NG))
    gas_count = 0
    for FILE in FILENAMES:

        WAVE, K_G = calc_k(FILE, VWAVE[0], VWAVE[-1], NLAYER, PRESS, TEMP)

        for IWAVE in range(len(WAVE)):
            for IG in range(NG):
                KOUT[IWAVE, :, gas_count, IG] = K_G[IWAVE, IG, :]

    return KOUT

def OVERLAP(delg, ng, k_g1, q1, k_g2, q2, g_ord):
    """
    Combine k distribution of two gases at a specified temperature and
    pressure and wavenumber.
    Inputs
    ------
    delg: array
        Widths of g-ordinate bins./or rather weights??
    ng: int
        Number of g-ordinates.
    k_g1: array
        k-distribution of first gas.
        k_g1 = kg_1(nwave, ng, npress, ntemp)
    q1: real
        Fraction of first gas.
    k_g2: array
        k-distribution of second gas.
    q2: real
        Fraction of second gas.
    g_ord: array
        The common g_ordinates.
    """
    # First deal with the case where one of the gas does not add opacity
    if k_g1[-1] == 0 or q1 == 0:
        k_g = k_g2 * q2/(q1+q2)
        return k_g
    if k_g2[-1] == 0 or q2 == 0:
        k_g = k_g1 * q1/(q1+q2)
        return k_g


    nloop = 0
    weight = np.zeros((ng*ng))
    contri =  np.zeros((ng*ng))
    for I in range(ng):
        for J in range(ng):
            weight[nloop] = delg[I]*delg[J]
            contri[nloop] = (k_g1[I]*q1 + k_g2[J]*q2)/(q1+q2)
            nloop+=1

    k_g = RANK(delg, ng, contri, weight, nloop, g_ord)

    return k_g

def RANK(delg, ng, contri, weight, nloop, g_ord):

    # form new g(k) by summing over weight. The new k(g) is then
    # found by ascending the arranged ks and getting the weighted averages
    # of the k values in the relevant g interval. Linear interpolation is
    # used to deal with steps crossing over the required g intervals.

    ascending_index = np.argsort(contri)
    contri = contri[ascending_index]
    weight = weight[ascending_index]
    gdist = np.concatenate([[0], np.cumsum(weight)])
    # gdist = np.cumsum(weight)
    """
    gdist = np.zeros(nloop+1)
    gdist[0] = 0
    gdist[1] = weight[0]
    for I in range(2,nloop):
        gdist[I] = weight[I] + gdist[I-1]
    """
    k_g = np.zeros(ng)

    ig = 0
    sum = 0.
    for I in range(nloop):
        # gdist holds cumulative distribution of ranked k
        if gdist[I] < g_ord[ig+1] and ig <= ng:
            k_g[ig] += contri[I]*weight[I]
            sum += weight[I]
        else:
            frac = (g_ord[ig+1]-gdist[I-1])/(gdist[I]-gdist[I-1])
            k_g[ig] += frac*contri[I]*weight[I]
            sum += frac*weight[I]
            k_g[ig] *= 1/sum
            ig += 1
            if ig <= ng:
                sum = (1-frac)*weight[I]
                k_g[ig] += (1-frac)*contri[I]*weight[I]


        #k_g[ng] = k_g[ng]/sum

    return k_g


def NGASCON(V, ID, ISO, AMOUNT,PP,PRESS,TEMP):
    """
    Inputs
    ------
    V0: real,
        Wavenumber.
    ID: integer array,
        Gas identifiers.
    ISO: integer array,
        Gas isotope identifiers.
    AMOUNT: real array,
        Absorber amount
    PP: real array,
        Partial pressure
    PRESS:
        Total pressure
    TEMP:
        Real

    Outputs
    -------
    AvgCONTMP: real
        Gas continuum abosprtion at V0
    """
    # currently no continuum
    AvgCONTMP = 0

    return AvgCONTMP

def NCIACON(V, PRESS, TEMP, INORMAL, ID, ISO, AMOUNT, PP):
    """
    Inputs
    ------
    V0: real,
        Wavenumber.
    ID: integer array,
        Gas identifiers.
    ISO: integer array,
        Gas isotope identifiers.
    AMOUNT: real array,
        Absorber amount
    PP: real array,
        Partial pressure
    PRESS:
        Total pressure
    TEMP:
        Real

    Outputs
    -------
    ABSORB: real,
        Calculated Optical depth
    """
    qh2=0
    qhe=0
    qn2=0
    qch4=0
    qco2=0

    NGAS = len(ID)
    # some routine to determine the aforementioned gas mixing ratios
    ABSORB = 0
    pass
    return ABSORB

def planck_wave(IWAVE, v, T):
    """
      Inputs
      ------
      IWAVE: int
          Indicates wavenumbers (0) or wavelength(1) for units.

      v: real
          Wavenumber (cm-1) or wavelength (um)

      T: real
          Temperature in K

      Returns
      -------
      Calculated function is in units of W cm-2 sr-1 cm for IWAVE=0
          or W cm-2 sr-1 um-1 for IWAVE=1
    """
    c1 = 1.1911e-12
    c2 = 1.439
    if IWAVE == 0:
        y = v
        a = c1*v**3
    else:
        y = 1e4/v
        a = c1*y**5/1e4
    if y == 0.0:
        planck_wave = 0.0
    else:
        tmp = c2*y/T
        b = np.exp(tmp) - 1
        planck_wave = a/b
    return planck_wave

import matplotlib.pyplot as plt

def find_nearest(array, value):
    """
    Find the closest value in an array
    INPUTS :
        array :: List of numbers
        value :: Value to search for
    OUTPUTS :
        closest_value :: Closest number to value in array
        index :: Index of closest_value within array
    CALLING SEQUENCE:
        closest_value,index = find_nearest(array,value)
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def read_ktahead(filename):
    """
    Read the header information in a correlated-k look-up
        table written with the standard format of Nemesis
    INPUTS :
        filename :: Name of the file (supposed to have a .kta extension)
    OUTPUTS :
        nwave :: Number of wavelength points
        vmin :: Minimum wavelength
        delv :: Spectral sampling
        npress :: Number of pressure levels
        ntemp :: Number of temperature levels
        gasID :: RADTRAN gas ID
        isoID :: RADTRAN isotopologue ID
        pressleves(np) :: Pressure levels (atm)
        templeves(np) :: Temperature levels (K)
    CALLING SEQUENCE:
        nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,\
            presslevels,templevels = read_ktahead(filename)
    """
    strlen = len(filename) # Opening file
    if filename[strlen-3:strlen] == 'kta':
        f = open(filename,'r')
    else:
        f = open(filename+'.kta','r')

    irec0 = np.fromfile(f,dtype='int32',count=1)
    nwave = np.fromfile(f,dtype='int32',count=1)
    vmin = np.fromfile(f,dtype='float32',count=1)
    delv = np.fromfile(f,dtype='float32',count=1)
    fwhm = np.fromfile(f,dtype='float32',count=1)
    npress = int(np.fromfile(f,dtype='int32',count=1))
    ntemp = int(np.fromfile(f,dtype='int32',count=1))
    ng = int(np.fromfile(f,dtype='int32',count=1))
    gasID = int(np.fromfile(f,dtype='int32',count=1))
    isoID = int(np.fromfile(f,dtype='int32',count=1))
    g_ord = np.fromfile(f,dtype='float32',count=ng)
    del_g = np.fromfile(f,dtype='float32',count=ng)
    presslevels = np.fromfile(f,dtype='float32',count=npress)
    N1 = abs(ntemp)
    if ntemp < 0:
        templevels = np.zeros([npress,N1])
        for i in range(npress):
            for j in range(N1):
                templevels[i,j] =  np.fromfile(f,dtype='float32',count=1)
    else:
        templevels = np.fromfile(f,dtype='float32',count=ntemp)

    return nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,\
        presslevels,templevels

def read_ktable(filename,wavemin,wavemax):
    """
    Read the correlated-k look-up table written in Nemesis format
    INPUTS :
        filename :: Name of the file (supposed to have a .kta extension)
        wavemin :: Minimum wavenumber to read (cm-1)
        wavemax :: Maximum wavenumber to read (cm-1)
    OUTPUTS :
        gasID :: Nemesis gas identifier
        isoID :: Nemesis isotopologue identifier
        nwave :: Number of wavenumbers
        wave(nwave) :: Wavenumbers or wavelengths
        fwhm :: Full width at half maximum
        ng :: Number of g-ordinates
        g_ord(ng) :: g-ordinates
        del_g(ng) :: Intervals of g-ordinates
        npress :: Number of pressure levels
        presslevels(npress) :: Pressure levels (atm)
        ntemp :: Number of temperature levels
        templevels(ntemp) :: Temperature levels (K)
        k_g(nwave,ng,npress,ntemp) :: K coefficients
    CALLING SEQUENCE:
        gasID,isoID,nwave,wave,fwhm,ng,g_ord,del_g,npress,presslevels,ntemp,\
            templevels,k_g = read_ktable(filename,wavemin,wavemax)
    """

    #Opening file
    strlen = len(filename)
    if filename[strlen-3:strlen] == 'kta':
        f = open(filename,'rb')
    else:
        f = open(filename+'.kta','rb')

    nbytes_int32 = 4
    nbytes_float32 = 4
    ioff = 0
    #Reading header
    irec0 = int(np.fromfile(f,dtype='int32',count=1))
    nwavekta = int(np.fromfile(f,dtype='int32',count=1))
    vmin = float(np.fromfile(f,dtype='float32',count=1))
    delv = float(np.fromfile(f,dtype='float32',count=1))
    fwhm = float(np.fromfile(f,dtype='float32',count=1))
    npress = int(np.fromfile(f,dtype='int32',count=1))
    ntemp = int(np.fromfile(f,dtype='int32',count=1))
    ng = int(np.fromfile(f,dtype='int32',count=1))
    gasID = int(np.fromfile(f,dtype='int32',count=1))
    isoID = int(np.fromfile(f,dtype='int32',count=1))

    ioff = ioff + 10 * nbytes_int32

    g_ord = np.zeros([ng])
    del_g = np.zeros([ng])
    templevels = np.zeros([ntemp])
    presslevels = np.zeros([npress])
    g_ord[:] = np.fromfile(f,dtype='float32',count=ng)
    del_g[:] = np.fromfile(f,dtype='float32',count=ng)

    ioff = ioff + 2*ng*nbytes_float32

    dummy = np.fromfile(f,dtype='float32',count=1)
    dummy = np.fromfile(f,dtype='float32',count=1)

    ioff = ioff + 2*nbytes_float32

    presslevels[:] = np.fromfile(f,dtype='float32',count=npress)
    templevels[:] = np.fromfile(f,dtype='float32',count=ntemp)

    ioff = ioff + npress*nbytes_float32+ntemp*nbytes_float32

    dummy = np.fromfile(f,dtype='float32',count=1)
    dummy = np.fromfile(f,dtype='float32',count=1)

    ioff = ioff + 2*nbytes_float32

    #Reading central wavelengths in non-uniform grid
    if delv>0.0:
        vmax = delv*nwavekta + vmin
        wavetot = np.linspace(vmin,vmax,nwavekta)
    else:
        wavetot = np.zeros([nwavekta])
        wavetot[:] = np.fromfile(f,dtype='float32',count=nwavekta)
        ioff = ioff + nwavekta*nbytes_float32

    #Calculating the wavenumbers to be read
    ins1 = np.where( (wavetot>=wavemin) & (wavetot<=wavemax) )
    ins = ins1[0]
    nwave = len(ins)
    wave = np.zeros([nwave])
    wave[:] = wavetot[ins]

    #Reading the k-coefficients
    k_g = np.zeros([nwave,ng,npress,ntemp])

    #Jumping until we get to the minimum wavenumber
    njump = npress*ntemp*ng*ins[0]
    ioff = njump*nbytes_float32 + (irec0-1)*nbytes_float32
    f.seek(ioff,0)

    #Reading the coefficients we require
    k_out = np.fromfile(f,dtype='float32',count=ntemp*npress*ng*nwave)
    il = 0
    for ik in range(nwave):
        for i in range(npress):
            for j in range(ntemp):
                k_g[ik,:,i,j] = k_out[il:il+ng]
                il = il + ng
    f.close()

    return gasID,isoID,nwave,wave,fwhm,ng,g_ord,del_g,\
        npress,presslevels,ntemp,templevels,k_g

def calc_k(filename,wavemin,wavemax,npoints,press,temp,MakePlot=False):

    """
    Calculate the k coefficients of a gas at a given pressure and
        temperature looking at pre-tabulated correlated-k tables
    INPUTS :
        filename :: Name of the file (supposed to have a .lta extension)
        nwave :: Number of wavenumbers (cm-1)
        wave :: Wavenumber (cm-1)
        npoints :: Number of p-T levels at which the absorption coefficient
            must be computed
        press(npoints) :: Pressure (atm)
        temp(npoints) :: Temperature (K)
    OPTIONAL INPUTS:
        MakePlot :: If True, a summary plot is generated
    OUTPUTS :
        wavek :: Calculation wavenumbers (cm-1)
        k(nwave,ng,npoints) :: K coefficients
    CALLING SEQUENCE:
        wavek,k = calc_k(filename,wavemin,wavemax,npoints,press,temp)
    """

    gasID,isoID,nwave,wave,fwhm,ng,g_ord,del_g,npress,presslevels,ntemp,\
        templevels,k_g = read_ktable(filename,wavemin,wavemax)

    #Interpolating to the correct pressure and temperature
    k_good = np.zeros([nwave,ng,npoints])
    for ipoint in range(npoints):
        press1 = press[ipoint]
        temp1 = temp[ipoint]

        #Getting the levels just above and below the desired points
        lpress  = np.log(press1)
        press0,ip = find_nearest(presslevels,press1)

        if presslevels[ip]>=press1:
            iphi = ip
            if ip==0:
                ipl = 0
            else:
                ipl = ip - 1
        elif presslevels[ip]<press1:
            ipl = ip
            if ip==npress-1:
                iphi = npress - 1
            else:
                iphi = ip + 1

        temp0,it = find_nearest(templevels,temp1)

        if templevels[it]>=temp1:
            ithi = it
            if it==0:
                itl = 0
            else:
                itl = it - 1
        elif templevels[it]<temp1:
            itl = it
            if it==ntemp-1:
                ithi = ntemp - 1
            else:
                ithi = it + 1

        plo = np.log(presslevels[ipl])
        phi = np.log(presslevels[iphi])
        tlo = templevels[itl]
        thi = templevels[ithi]
        klo1 = np.zeros([nwave,ng])
        klo2 = np.zeros([nwave,ng])
        khi1 = np.zeros([nwave,ng])
        khi2 = np.zeros([nwave,ng])
        klo1[:] = k_g[:,:,ipl,itl]
        klo2[:] = k_g[:,:,ipl,ithi]
        khi2[:] = k_g[:,:,iphi,ithi]
        khi1[:] = k_g[:,:,iphi,itl]

        #Interpolating to get the k-coefficients at desired p-T
        if ipl==iphi:
            v = 0.5
        else:
            v = (lpress-plo)/(phi-plo)

        if itl==ithi:
            u = 0.5
        else:
            u = (temp1-tlo)/(thi-tlo)

        k_good[:,:,ipoint] = (1.0-v)*(1.0-u)*klo1[:,:] \
            + v*(1.0-u)*khi1[:,:] + v*u*khi2[:,:] + (1.0-v)*u*klo2[:,:]

    """
    if MakePlot==True:
        fig, ax = plt.subplots(1,1,figsize=(10,6))
        k_abs = np.matmul(k_good[:,:,npoints-1], del_g)
        k_abslo1 = np.matmul(klo1[:,:], del_g)
        k_abslo2 = np.matmul(klo2[:,:], del_g)
        k_abshi1 = np.matmul(khi1[:,:], del_g)
        k_abshi2 = np.matmul(khi2[:,:], del_g)
        ax.semilogy(wave,k_abslo1,label='p = '+str(np.exp(plo))+' atm - T = '+str(tlo)+' K')
        ax.semilogy(wave,k_abslo2,label='p = '+str(np.exp(plo))+' atm - T = '+str(thi)+' K')
        ax.semilogy(wave,k_abshi1,label='p = '+str(np.exp(phi))+' atm - T = '+str(tlo)+' K')
        ax.semilogy(wave,k_abshi2,label='p = '+str(np.exp(phi))+' atm - T = '+str(thi)+' K')
        ax.semilogy(wave,k_abs,label='p = '+str(press1)+' atm - T = '+str(temp1)+' K',color='black')
        ax.legend()
        ax.grid()
        plt.tight_layout()
        plt.show()
    """
    return wave,k_good

wasp = np.array([3.10545e+25,
                 2.99568e+25,
                 2.82785e+25,
                 2.73007e+25,
                 2.69301e+25,
                 2.59410e+25,
                 2.47917e+25,
                 2.45594e+25,
                 2.31012e+25,
                 2.34431e+25,
                 2.29190e+25,
                 2.20532e+25,
                 2.15548e+25,
                 1.23401e+24,
                 4.42220e+23])



H = np.array([0.,178.74 ,333.773,460.83 ,572.974,680.655,
        787.549,894.526,1001.764,1109.3  ])*1e3

P = np.array([1.9739e+01,3.9373e+00,7.8539e-01,1.5666e-01,3.1250e-02,
       6.2336e-03,1.2434e-03,2.4804e-04,4.9476e-05,9.8692e-06])*101325

T = np.array([1529.667,1408.619,1128.838,942.708,879.659,864.962,
        861.943,861.342,861.222,861.199])


# H2O,CO2,He,H2
RADIUS = 74065.70e3
ID, ISO = [1,2,40,39], [0,0,0,0]
NVMR = len(ISO)
VMR = np.array([[0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915],
                [0.001,0.14985,0.84915]])

PP = np.array([[1.38394260e+03, 2.07383799e+05, 1.17517486e+06],
               [2.90141092e+02, 4.34776426e+04, 2.46373308e+05],
               [6.38644910e+01, 9.57009398e+03, 5.42305325e+04],
               [1.44413798e+01, 2.16404076e+03, 1.22628976e+04],
               [3.34426797e+00, 5.01138555e+02, 2.83978515e+03],
               [7.99615911e-01, 1.19822444e+02, 6.78993851e+02],
               [1.95720114e-01, 2.93286591e+01, 1.66195735e+02],
               [4.84525389e-02, 7.26061295e+00, 4.11434734e+01],
               [1.20207195e-02, 1.80130481e+00, 1.02073939e+01],
               [2.97111025e-03, 4.45220871e-01, 2.52291827e+00]])

AMOUNT = np.array([[1.03055120e+28, 1.54428097e+30, 8.75092549e+30],
                    [2.15916562e+27, 3.23550968e+29, 1.83345548e+30],
                    [4.73939916e+26, 7.10198964e+28, 4.02446080e+29],
                    [1.05719726e+26, 1.58421009e+28, 8.97719053e+28],
                    [2.36105206e+25, 3.53803652e+27, 2.00488736e+28],
                    [5.40063511e+24, 8.09285172e+26, 4.58594931e+27],
                    [1.27266617e+24, 1.90709026e+26, 1.08068448e+27],
                    [3.08166897e+23, 4.61788095e+25, 2.61679921e+26],
                    [7.64558085e+22, 1.14569029e+25, 6.49224498e+25],
                    [1.93872185e+22, 2.90517470e+24, 1.64626566e+25]])

SCALE = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

NTOTAM = np.array([1.0312e+27, 2.1605e+26, 4.7421e+25, 1.0578e+25, 2.3625e+24,
       5.4042e+23, 1.2735e+23, 3.0836e+22, 7.6502e+21, 1.9399e+21])

NPRESS = np.array([1.38394260e+06, 2.90141092e+05, 6.38644910e+04, 1.44413798e+04,
       3.34426797e+03, 7.99615911e+02, 1.95720114e+02, 4.84525389e+01,
       1.20207195e+01, 2.97111025e+00])

NTEMP = np.array([1483.088, 1311.001, 1076.866,  932.43 ,  880.159,  865.926,
        862.359,  861.481,  861.264,  861.211])

NTEMP = NTEMP

VWAVE =  np.array([1., 1.24749994, 1.28250003, 1.3175    , 1.35249996,
       1.38750005, 1.42250001, 1.45749998, 1.49249995, 1.52750003,
       1.5625    , 1.59749997, 1.63250005, 3.5999999 , 5       ])

DELG = np.array([0.008807  , 0.02030071, 0.03133602, 0.04163837, 0.05096506,
       0.05909727, 0.06584432, 0.07104805, 0.0745865 , 0.07637669,
       0.07637669, 0.0745865 , 0.07104805, 0.06584432, 0.05909727,
       0.05096506, 0.04163837, 0.03133602, 0.02030071, 0.008807  ])

G_ORD = np.array([0.0034357 , 0.01801404, 0.04388279, 0.08044151, 0.12683405,
       0.18197316, 0.2445665 , 0.31314695, 0.38610709, 0.46173674,
       0.53826326, 0.61389291, 0.68685305, 0.7554335 , 0.81802684,
       0.87316597, 0.91955847, 0.95611721, 0.98198599, 0.99656433])

FILENAMES = ['h2o']

"""
R = IMOD3(RADIUS=RADIUS, ID=ID, ISO=ISO, TOTAM=NTOTAM, PRESS=NPRESS, TEMP=NTEMP,
        PP=PP, AMOUNT=AMOUNT, SCALE=SCALE,
        VWAVE=VWAVE, FILENAMES=FILENAMES, DELG=DELG, G_ORD=G_ORD)
print(R)
"""
RADIUS=RADIUS
ID=ID
ISO=ISO
TOTAM=NTOTAM
PRESS=NPRESS
TEMP=NTEMP*2
PP=PP
AMOUNT=AMOUNT
SCALE=SCALE
VWAVE=VWAVE
FILENAMES=FILENAMES
DELG=DELG
G_ORD=G_ORD



FRAC = (PP.T/PRESS).T
print('FRAC', FRAC)
UTOTL = TOTAM*1e27
print(UTOTL)
NG = len(G_ORD)
print('NG',NG)
NWAVE = len(VWAVE)
print('NWAVE', NWAVE)
NLAYER = len(PRESS)
print('PRESS', PRESS)
NGAS = len(FILENAMES)
print('NGAS',NGAS)

# Interpolate k-tables for k values at each wave bin, for each layer,
# from each gas and at each g ordinate.
# K_DATA.shape = (NWAVE, NLAYER, NGAS, NG)
K_DATA = GET_K(FILENAMES,NG,VWAVE,PRESS,TEMP)

# initialise output array
RADIANCE = np.zeros(NWAVE)

for I in range(NWAVE):
    print('IWAVE', I+1)
    """
    if ISPACE == 1:
        # wavenumber space
        VV = 1e4/VWAVE[I]
    else:
        # wavelength space
        VV = VWAVE[I]
    """
    VV = VWAVE[I]

    # Get k data for the i-th layer
    K_DATA_I = K_DATA[I]
    # combined k data from different gas assuming random overlaps
    K_G = np.zeros(NG)
    K_G2 = np.zeros(NG)
    KL_G = np.zeros((NG, NLAYER))

    for J in range(NLAYER):

        for K in range(NGAS):
            K_G2 = K_DATA_I[J,K,:]

            if K == 0:
                q1 = FRAC[J,K]
                #print('q1',q1)
                K_G = K_G2
                #print(K_G)

            else:
                q2 = FRAC[J,K]
                K_G1 = K_G
                if q2>0 and q1>0:
                    K_G = OVERLAP(DELG, NG, K_G1, q1, K_G2, q2, G_ORD)
                else:
                    if q2>0:
                        K_G = K_G2
                q1 = q1 + q2
        print(K_G)
        KL_G[:,J] = K_G
        print(KL_G)
        
        # print('KL_G', KL_G)
    CORKOUT = np.zeros(NG)
    TAUCON = np.zeros(NLAYER)
    for IG in range(NG):
        print("IG", IG)
        TAUTMP = TAUCON + KL_G[IG,:] * UTOTL
        TAUS = TAUTMP * SCALE
        # print('TAUS', TAUS)
        XFAC = 1

    

        # blackbody radiation for each layer
        BB = planck_wave(IWAVE=1,v=VV,T=TEMP)
        # print(BB)
        # TR: Cumulative transmission function
        # Radiance = sum(Blackbody radiation * Transmission Intervals)
        TR = np.exp(-np.cumsum(TAUS))
        TR = np.concatenate([[1],TR])
        DEL_TR = -(TR[1:]-TR[:-1])
        CORKOUT[IG] = np.sum(BB * DEL_TR)
        print(CORKOUT[IG])

    # Integrate over g ordinates for each wave bin
    # len(CORKOUT) = len(DELG) = NG
    RADIANCE[I] = np.sum(CORKOUT * DELG)
    
o = RADIANCE/wasp
plt.plot(VWAVE,o)