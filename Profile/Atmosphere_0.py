#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:27:12 2021

@author: jingxuanyang

Basic Atmosphere Class.
"""
import numpy as np

class Atmosphere_0:
    """
    Clear atmosphere. Simplest possible profile.
    """
    def __init__(self, runname='wasp43b', NP=10, NVMR=6, ID=[0,0,0,0,0,0],
                ISO=[0,0,0,0,0,0], LATITUDE=0.0, IPLANET=1, AMFORM=1):
        """
        Set up an atmosphere profile with NP points and NVMR gases.
        Use the class methods to edit Height, Pressure, Temperature and
        gas Volume Mixing Ratios after creating an Atmosphere_0 object, e.g.
            atm = Atmosphere_0('Jupiter',100,6,[1,2,4,5,39,40],
                      [0,0,0,0,0,0],0,1,1)
            atm.edit_H(your_height_profile)
            atm.editP(your_pressure_profile)
            atm.edit_T(your_temperature_profile)
            atm.edit_VMR(your_VMR_profile)
        Can write to a *.prf file in Nemesis format.

        Inputs
        ------
        @param runname: str
            Name of this particular profile (no space)
        @param NP: int,
            Number of points defined in the profile
        @param NVMR: int
            Number of gases in the profile (expressed as volume mixing ratios)
        @param ID: 1D array
            Gas ID for each gas to be defined in the profile
        @param ISO: 1D array
            Isotope ID for each gas, default 0 for all
            isotopes in terrestrial relative abundance
        @param IPLANET: int
            Planet ID
        @param LATITUDE: real
            Planetocentric latitude
        @param AMFORM: int,
            Format of the atmospheric profile, default AMFORM=1:
            assume that at each level the VMRs add up to 1.0.

        Attributes
        ----------
        @attribute H: 1D array
            Height in m of each points above reference level:
            ground for terrestrial planets,
            usually 1-bar level for gas giants
        @attribute P: 1D array
            Pressure in Pa of each point defined in the profile
        @attribute T: 1D array
            Temperature in K at each point defined in the profile
        @attribute VMR: 2D array
            VMR[i,j] is Volume Mixing Ratio of gas j at vertical point i
            the column j corresponds to the gas with RADTRANS ID ID[j]
            and RADTRANS isotope ID ISO[j]

        Methods
        -------
        Atmosphere_0.edit_H
        Atmosphere_0.edit_P
        Atmosphere_0.edit_T
        Atmosphere_0.edit_VMR
        Atmosphere_0.write_to_file
        Atmosphere_0.check
        """

        assert type(runname) == str and len(runname) <= 100,\
            'runname should be a string <= 100 char'
        assert type(NP) == int and type(NVMR) == int,\
            'NP and NVMR should be integers'
        try:
            assert type(ID) == np.ndarray
        except:
            assert type(ID) == list, 'ID should be an array'
            ID = np.array(ID)
        try:
            assert type(ISO) == np.ndarray
        except:
            assert type(ISO) == list, 'ISO should be an array'
            ISO = np.array(ISO)
        assert len(ID) == NVMR, "len(ID) should be equal to NVMR"
        for i in range(NVMR):
            assert type(ID[i]) == np.int64,\
                'RADTRANS gas IDs should be integers'
            assert type(ISO[i]) == np.int64,\
                'RADTRANS isotope IDs should be integers'
        assert len(ID) == len(ISO),\
            "len(ID) should be equal to len(ISO)"
        assert type(IPLANET) == int and type(AMFORM) == int,\
            'IPLANET and AMFORM should be integers'
        assert abs(LATITUDE) <= 90, 'LATITUDE should be <= 90 deg'

        self.runname = runname
        self.NP = NP
        self.NVMR = NVMR
        self.ID = ID
        self.ISO = ISO
        self.IPLANET = IPLANET
        self.LATITUDE = LATITUDE
        self.AMFORM = AMFORM

        # Input the following profiles using the edit_ methods.
        self.H = None # np.zeros(NP)
        self.P = None # np.zeros(NP)
        self.T =  None # np.zeros(NP)
        self.VMR = None # np.zeros((NP, NVMR))

    def edit_H(self, H_array):
        """
        Edit the Height profile.
        @param H_array: 1D array
            Heights of the vertical points in m
        """
        H_array = np.array(H_array)
        assert len(H_array) == self.NP, 'H should have NP elements'
        assert ((H_array[1:]-H_array[:-1])>0).all(),\
            'H should be strictly increasing'
        self.H = H_array

    def edit_P(self, P_array):
        """
        Edit the Pressure profile.
        @param P_array: 1D array
            Pressures of the vertical points in Pa
        """
        P_array = np.array(P_array)
        assert len(P_array) == self.NP, 'P should have NP elements'
        assert (P_array>0).all() == True, 'P should be positive'
        assert ((P_array[1:]-P_array[:-1])<0).all(),\
            'P should be strictly decreasing'
        self.P = P_array

    def edit_T(self, T_array):
        """
        Edit the Temperature profile.
        @param T_array: 1D array
            Temperature of the vertical points in K
        """
        T_array = np.array(T_array)
        assert len(T_array) == self.NP, 'T should have NP elements'
        assert (T_array>0).all() == True, 'T should be positive'
        self.T = T_array

    def edit_VMR(self, VMR_array):
        """
        Edit the gas Volume Mixing Ratio profile.
        @param VMR_array: 2D array
            NP by NVMR array containing the Volume Mixing Ratio of gases.
            VMR_array[i,j] is the volume mixing ratio of gas j at point i.
        """
        VMR_array = np.array(VMR_array)
        try:
            assert VMR_array.shape == (self.NP, self.NVMR),\
                'VMR should be NP by NVMR.'
        except:
            assert VMR_array.shape == (self.NP,) and self.NVMR==1,\
                'VMR should be NP by NVMR.'
        assert (VMR_array>=0).all() == True,\
            'VMR should be non-negative'
        self.VMR = VMR_array

    def write_to_file(self):
        """
        Write the current profile to a runname.prf file in Nemesis format.
        """
        self.check()
        f = open('{}.prf'.format(self.runname),'w')
        f.write('{}\n'.format(self.AMFORM))
        f.write('{:<10} {:<10} {:<10} {:<10}\n'
                .format(self.IPLANET, self.LATITUDE, self.NP, self.NVMR))
        for i in range(self.NVMR):
            f.write('{:<5} {:<5}\n'.format(self.ID[i], self.ISO[i]))
        f.write('{:<15} {:<15} {:<15} '
                .format('H(km)', 'press(atm)', 'temp(K)'))
        for i in range(self.NVMR):
            f.write('{:<15} '.format('VMR gas{}'.format(i+1)))
        for i in range(self.NP):
            f.write('\n{:<15.3f} {:<15.5E} {:<15.3f} '
                    .format(self.H[i]*1e-3, self.P[i]/101325, self.T[i]))
            for j in range(self.NVMR):
                f.write('{:<15.5E} '.format(self.VMR[i][j]))
        f.close()

    def check(self):
        assert isinstance(self.H, np.ndarray),\
            'Need to input height profile'
        assert isinstance(self.P, np.ndarray),\
            'Need to input pressure profile'
        assert isinstance(self.T, np.ndarray),\
            'Need to input temperature profile'
        assert isinstance(self.VMR, np.ndarray),\
            'Need to input VMR profile'
        return True


atm0 = Atmosphere_0()
for i in range(1):
    #atm0 = Atmosphere_0()
    #atm0.write_to_file()

    # create profiles from external models
    H = np.linspace(0,9000,10)
    P = np.logspace(6,1,10)
    T = np.linspace(40,20,10)**2
    VMR = np.array([np.ones(10)*1.6e-6,
                          np.ones(10)*1.6e-6,
                          np.ones(10)*1.6e-6,
                          np.ones(10)*1.6e-6,
                          np.ones(10)*1.6e-6,
                          np.ones(10)*1.6e-6,]).T

    # add profiles to atmosphere
    atm0.edit_H(H)
    atm0.edit_P(P)
    atm0.edit_T(T)
    atm0.edit_VMR(VMR)
    atm0.check()
    #atm0.write_to_file()
