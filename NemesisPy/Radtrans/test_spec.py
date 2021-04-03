# -*- coding: utf-8 -*-

from spectroscopy import read_ktahead, read_ktable, calc_k
import numpy as np
# k(nwave,ng,npoints)
a, b = calc_k('h2o',1.0,4.5,3,np.array([1,0.1,0.01]),np.array([2000,1000,500]),MakePlot=False)

nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,presslevels,templevels \
    = read_ktahead('co2')

#def calc_k(filename,wavemin,wavemax,npoints,press,temp,MakePlot=False):
wave,k = calc_k('co2',1,5,1,[1],[1000])

gasID,isoID,nwave,wave,fwhm,ng,g_ord,del_g,npress,presslevels,ntemp,templevels,k_g = read_ktable('co2',1,5)


def rank(g_ord, weight, k_g):
    ng = len(g_ord)

    
"""
Parse 'co2.kta' 
nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,presslevels,templevels \
    = read_ktahead('co2')
    nwave = array([17], dtype=int32)
    vmin = array([1.1425], dtype=float32)
    delv = array([-1.], dtype=float32)
    fwhm = array([0.], dtype=float32)
    npress = 20
    ntemp = 20
    ng = 20
    gasID = 1
    isoID = 0
    
    g_ord = array([0.0034357 , 0.01801404, 0.04388279, 0.08044151, 0.12683405,
       0.18197316, 0.2445665 , 0.31314695, 0.3861071 , 0.46173674,
       0.53826326, 0.6138929 , 0.68685305, 0.7554335 , 0.81802684,
       0.87316597, 0.91955847, 0.9561172 , 0.981986  , 0.9965643 ],
      dtype=float32)
    
    del_g = array([0.008807  , 0.02030071, 0.03133602, 0.04163837, 0.05096506,
       0.05909727, 0.06584432, 0.07104805, 0.0745865 , 0.07637669,
       0.07637669, 0.0745865 , 0.07104805, 0.06584432, 0.05909727,
       0.05096506, 0.04163837, 0.03133602, 0.02030071, 0.008807  ],
      dtype=float32)
    
    presslevels = array([0.00000000e+00, 0.00000000e+00, 3.05902319e-07, 8.58439876e-07,
       2.40900340e-06, 6.76027776e-06, 1.89710809e-05, 5.32376871e-05,
       1.49398518e-04, 4.19250515e-04, 1.17652444e-03, 3.30162724e-03,
       9.26521700e-03, 2.60005593e-02, 7.29641989e-02, 2.04756334e-01,
       5.74598432e-01, 1.61247122e+00, 4.52500534e+00, 1.26983185e+01],
      dtype=float32)
    
    templevels = array([  35.63472,  100.00029,  100.     ,  250.     ,  400.     ,
        550.     ,  700.     ,  850.     , 1000.     , 1150.     ,
       1300.     , 1450.     , 1600.     , 1750.     , 1900.     ,
       2050.     , 2200.     , 2350.     , 2500.     , 2650.     ],
      dtype=float32)
"""
