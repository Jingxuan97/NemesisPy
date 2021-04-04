# -*- coding: utf-8 -*-

from spectroscopy import read_ktahead, read_ktable, calc_k
import numpy as np

"""
nwave,vmin,delv,fwhm,npress,ntemp,ng,gasID,isoID,g_ord,del_g,\
    presslevels(npk),templevels(ntk) = read_ktahead(filename)
    
gasID,isoID,nwave,wave,fwhm,ng,g_ord,del_g,npress,presslevels,ntemp,templevels,\
    k_g(nwave,ng,npress,ntemp)= read_ktable(filename,wavemin,wavemax)

wavek,k = calc_k(filename,wavemin,wavemax,npoints,press,temp)
"""