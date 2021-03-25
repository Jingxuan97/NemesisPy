#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
from Profile import Atmosphere_1
from Layer import Layer
import matplotlib.pyplot as plt
###############################################################################
#                                                                             #
#                               MODEL input                                   #
#                                                                             #
###############################################################################
H = np.array([0.,178.74 ,333.773,460.83 ,572.974,680.655,
        787.549,894.526,1001.764,1109.3  ])*1e3

P = np.array([1.9739e+01,3.9373e+00,7.8539e-01,1.5666e-01,3.1250e-02,
       6.2336e-03,1.2434e-03,2.4804e-04,4.9476e-05,9.8692e-06])*101325

T = np.array([1529.667,1408.619,1128.838,942.708,879.659,864.962,
        861.943,861.342,861.222,861.199])

# H2O,He,H2
ID, ISO = [1,40,39], [0,0,0]
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

NP  = 10
runname, LATITUDE, IPLANET, AMFORM = 'wasp43b', 0.0, 0, 1
Atm = Atmosphere_1(runname=runname,NP=NP,NVMR=NVMR,ID=ID,ISO=ISO,
                    LATITUDE=LATITUDE,IPLANET=IPLANET,AMFORM=AMFORM)
Atm.edit_H(H)
Atm.edit_P(P)
Atm.edit_T(T)
Atm.edit_VMR(VMR)
###############################################################################
#                                                                             #
#                               Layer input                                   #
#                                                                             #
###############################################################################
# Calculate average layer properties
NLAY = 10
RADIUS = 74065.70e3
LAYANG = 0
LAYINT = 1
H,P,T = Atm.H, Atm.P, Atm.T
VMR = Atm.VMR
LAYTYP, AMFORM, LAYHT = 1, 1, 0.0
H_base, P_base, INTERTYP = None, None, 1

# layer properties
Layer = Layer(RADIUS=RADIUS, LAYTYP=LAYTYP, NLAY=NLAY, LAYINT=LAYINT, NINT=101,
              AMFORM=AMFORM, INTERTYP=INTERTYP, H_base=H_base, P_base=P_base)

BASEH, BASEP, BASET, HEIGHT, PRESS, TEMP, TOTAM, AMOUNT, PP, LAYSF, DELH\
    = Layer.integrate(H=H,P=P,T=T, LAYANG=LAYANG, ID=ID,VMR=VMR)


print('{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} '.format(
'BASEH', 'DELH', "BASEP", 'BASET', 'TOTAM', 'PRESS', 'TEMP'),end="")
for j in range(NVMR):
    print('{:<8} {:<8}'.format('AMOUNT', 'PP'), end='')
for i in range(len(BASEH)):
    print('\n{} {:<8.2f} {:<8.2f} {:<8.3E} {:<8.0f} {:<8.3e} {:<8.3e} {:<8.0f} '.format(i+1,
    BASEH[i]/1e3, DELH[i]/1e3, BASEP[i]/101325, BASET[i],
    TOTAM[i]/1e4, PRESS[i]/101325, TEMP[i]),end="")
    for j in range(NVMR):
        print('{:<10.3E} {:<10.3E}'.format(AMOUNT[i,j]/1e4,PP[i,j]/101325),end="")
    print('\n')
###############################################################################
#                                                                             #
#                               Compare                                       #
#                                                                             #
###############################################################################
NBASEH = np.array([   0.  ,  170.95,  319.07,  441.14,  547.66,  647.41,  744.11,
        838.79,  931.45, 1021.76])

NBASEP = np.array([1.9739e+01, 4.6262e+00, 1.0842e+00, 2.5410e-01, 5.9554e-02,
       1.3957e-02, 3.2711e-03, 7.6665e-04, 1.7968e-04, 4.2110e-05])

NBASET = np.array([1529.667, 1413.896, 1155.364,  971.556,  893.888,  869.5  ,
        863.17 ,  861.655,  861.301,  861.218])

NTOTAM = np.array([1.0312e+27, 2.1605e+26, 4.7421e+25, 1.0578e+25, 2.3625e+24,
       5.4042e+23, 1.2735e+23, 3.0836e+22, 7.6502e+21, 1.9399e+21])

NPRESS = np.array([1.3659e+01, 2.8635e+00, 6.3027e-01, 1.4252e-01, 3.3006e-02,
       7.8920e-03, 1.9317e-03, 4.7819e-04, 1.1864e-04, 2.9322e-05])

NTEMP = np.array([1483.088, 1311.001, 1076.866,  932.43 ,  880.159,  865.926,
        862.359,  861.481,  861.264,  861.211])

plt.plot(TOTAM*1e-4, BASEH*1e-3, label='Python')
# plt.scatter(NTOTAM, NBASEH, label='Fortran', color='k')
plt.ylabel('height (km)')
plt.xlabel('total amount (no./cm^2)')
plt.legend()


###############################################################################
#                                                                             #
#                               Data                                          #
#                                                                             #
###############################################################################
"""NEMESIS
  1    0.00  170.95 0.19739E+02 1529.667 0.10312E+28 0.13659E+02 1483.088
    0.10312E+25 0.13659E-01 0.15452E+27 0.20467E+01 0.87563E+27 0.11598E+02

  2  170.95  148.13 0.46262E+01 1413.896 0.21605E+27 0.28635E+01 1311.001
    0.21605E+24 0.28635E-02 0.32375E+26 0.42910E+00 0.18346E+27 0.24315E+01

  3  319.07  122.06 0.10842E+01 1155.364 0.47421E+26 0.63027E+00 1076.866
    0.47421E+23 0.63027E-03 0.71061E+25 0.94445E-01 0.40268E+26 0.53519E+00

  4  441.14  106.53 0.25410E+00  971.556 0.10578E+26 0.14252E+00  932.430
    0.10578E+23 0.14252E-03 0.15851E+25 0.21357E-01 0.89825E+25 0.12102E+00

  5  547.66   99.74 0.59554E-01  893.888 0.23625E+25 0.33006E-01  880.159
    0.23625E+22 0.33006E-04 0.35403E+24 0.49459E-02 0.20061E+25 0.28027E-01

  6  647.41   96.70 0.13957E-01  869.500 0.54042E+24 0.78920E-02  865.926
    0.54042E+21 0.78920E-05 0.80982E+23 0.11826E-02 0.45890E+24 0.67015E-02

  7  744.11   94.68 0.32711E-02  863.170 0.12735E+24 0.19317E-02  862.359
    0.12735E+21 0.19317E-05 0.19083E+23 0.28947E-03 0.10814E+24 0.16403E-02

  8  838.79   92.66 0.76665E-03  861.655 0.30836E+23 0.47819E-03  861.481
    0.30836E+20 0.47819E-06 0.46207E+22 0.71656E-04 0.26184E+23 0.40605E-03

  9  931.45   90.32 0.17968E-03  861.301 0.76502E+22 0.11864E-03  861.264
    0.76502E+19 0.11864E-06 0.11464E+22 0.17778E-04 0.64962E+22 0.10074E-03

 10 1021.76   87.54 0.42110E-04  861.218 0.19399E+22 0.29322E-04  861.211
    0.19399E+19 0.29322E-07 0.29070E+21 0.43939E-05 0.16473E+22 0.24899E-04
"""

"""NEMESIS
 data = np.array([[1.000000e+00, 0.000000e+00, 1.709500e+02, 1.973900e+01,
        1.529667e+03, 1.031200e+27, 1.365900e+01, 1.483088e+03,
        1.031200e+24, 1.365900e-02, 1.545200e+26, 2.046700e+00,
        8.756300e+26, 1.159800e+01],
       [2.000000e+00, 1.709500e+02, 1.481300e+02, 4.626200e+00,
        1.413896e+03, 2.160500e+26, 2.863500e+00, 1.311001e+03,
        2.160500e+23, 2.863500e-03, 3.237500e+25, 4.291000e-01,
        1.834600e+26, 2.431500e+00],
       [3.000000e+00, 3.190700e+02, 1.220600e+02, 1.084200e+00,
        1.155364e+03, 4.742100e+25, 6.302700e-01, 1.076866e+03,
        4.742100e+22, 6.302700e-04, 7.106100e+24, 9.444500e-02,
        4.026800e+25, 5.351900e-01],
       [4.000000e+00, 4.411400e+02, 1.065300e+02, 2.541000e-01,
        9.715560e+02, 1.057800e+25, 1.425200e-01, 9.324300e+02,
        1.057800e+22, 1.425200e-04, 1.585100e+24, 2.135700e-02,
        8.982500e+24, 1.210200e-01],
       [5.000000e+00, 5.476600e+02, 9.974000e+01, 5.955400e-02,
        8.938880e+02, 2.362500e+24, 3.300600e-02, 8.801590e+02,
        2.362500e+21, 3.300600e-05, 3.540300e+23, 4.945900e-03,
        2.006100e+24, 2.802700e-02],
       [6.000000e+00, 6.474100e+02, 9.670000e+01, 1.395700e-02,
        8.695000e+02, 5.404200e+23, 7.892000e-03, 8.659260e+02,
        5.404200e+20, 7.892000e-06, 8.098200e+22, 1.182600e-03,
        4.589000e+23, 6.701500e-03],
       [7.000000e+00, 7.441100e+02, 9.468000e+01, 3.271100e-03,
        8.631700e+02, 1.273500e+23, 1.931700e-03, 8.623590e+02,
        1.273500e+20, 1.931700e-06, 1.908300e+22, 2.894700e-04,
        1.081400e+23, 1.640300e-03],
       [8.000000e+00, 8.387900e+02, 9.266000e+01, 7.666500e-04,
        8.616550e+02, 3.083600e+22, 4.781900e-04, 8.614810e+02,
        3.083600e+19, 4.781900e-07, 4.620700e+21, 7.165600e-05,
        2.618400e+22, 4.060500e-04],
       [9.000000e+00, 9.314500e+02, 9.032000e+01, 1.796800e-04,
        8.613010e+02, 7.650200e+21, 1.186400e-04, 8.612640e+02,
        7.650200e+18, 1.186400e-07, 1.146400e+21, 1.777800e-05,
        6.496200e+21, 1.007400e-04],
       [1.000000e+01, 1.021760e+03, 8.754000e+01, 4.211000e-05,
        8.612180e+02, 1.939900e+21, 2.932200e-05, 8.612110e+02,
        1.939900e+18, 2.932200e-08, 2.907000e+20, 4.393900e-06,
        1.647300e+21, 2.489900e-05]])
"""

"""
prf = np.array([[ 0.000,0.19739E+02,1529.6670,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 178.740,0.39373E+01,1408.6190,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 333.773,0.78539E+00,1128.8380,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 460.830,0.15666E+00,942.7080,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 572.974,0.31250E-01,879.6590,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 680.655,0.62336E-02,864.9620,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 787.549,0.12434E-02,861.9430,0.10000E-02,0.14985E+00,0.84915E+00],
      [ 894.526,0.24804E-03,861.3420,0.10000E-02,0.14985E+00,0.84915E+00],
      [1001.764,0.49476E-04,861.2220,0.10000E-02,0.14985E+00,0.84915E+00],
      [1109.300,0.98692E-05,861.1990,0.10000E-02,0.14985E+00,0.84915E+00],])
"""