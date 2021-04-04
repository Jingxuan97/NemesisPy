import numpy as np
radtran_ID = {
    'H2O': 1,
    'CO2': 2,
    'O3': 3,
    'N2O': 4,
    'CO': 5,
    'CH4': 6,
    'O2': 7,
    'NO': 8,
    'SO2': 9,
    'NO2': 10,
    'NH3': 11,
    'HNO3': 12,
    'OH': 13,
    'HF': 14,
    'HCL': 15,
    'HBr': 16,
    'HI': 17,
    'ClO': 18,
    'OCS': 19,
    'H2CO': 20,
    'HOCl': 21,
    'N2': 22,
    'HCN': 23,
    'CH3Cl': 24,
    'H2O2': 25,
    'C2H2': 26,
    'C2H6': 27,
    'PH3': 28,
    'C2N2': 29,
    'C4H2': 30,
    'HC3N': 31,
    'C2H4': 32,
    'GeH4': 33,
    'C3H8': 34,
    'HCOOH': 35,
    'H2S': 36,
    'COF2': 37,
    'SF6': 38,
    'H2': 39,
    'He': 40,
    'AsH3': 41,
    'C3H4': 42,
    'ClONO2': 43,
    'HO2': 44,
    'O': 45,
    'NO+': 46,
    'CH3OH': 47,
    'H': 48,
    'C6H6': 49,
    'CH3CN': 50,
    'CH2NH': 51,
    'C2H3CN': 52,
    'HCP': 53,
    'CS': 54,
    'HC5N': 55,
    'HC7N': 56,
    'C2H5CN': 57,
    'CH3NH2': 58,
    'HNC': 59,
    'Na': 60,
    'K': 61,
    'TiO': 62,
    'VO': 63,
    'CH2CCH2': 64,
    'C4N2': 65,
    'C5H5N': 66,
    'C5H4N2': 67,
    'C7H8': 68,
    'C8H6': 69,
    'C5H5CN': 70,
    'HOBr': 71,
    'CH3Br': 72,
    'CF4': 73,
    'SO3': 74,
    'Ne': 75,
    'Ar': 76,
    'COCl2': 77,
    'SO': 78,
    'H2SO4': 79,
    'e–': 80,
    'H3+': 81,
    'FeH': 82,
    'AlO': 83,
    'AlCl': 84,
    'AlF': 85,
    'AlH': 86,
    'BeH': 87,
    'C2': 88,
    'CaF': 89,
    'CaH': 90,
    'H–': 91,
    'CaO': 92,
    'CH': 93,
    'CH3': 94,
    'CH3F': 95,
    'CN': 96,
    'CP': 97,
    'CrH': 98,
    'HD+': 99,
    'HeH+': 100,
    'KCl': 101,
    'KF': 102,
    'LiCl': 103,
    'LiF': 104,
    'LiH': 105,
    'LiH+': 106,
    'MgF': 107,
    'MgH': 108,
    'MgO': 109,
    'NaCl': 110,
    'NaF': 111,
    'NaH': 112,
    'NH': 113,
    'NS': 114,
    'OH+': 115,
    'cis-P2H2': 116,
    'trans-P2H2': 117,
    'PH': 118,
    'PN': 119,
    'PO': 120,
    'PS': 121,
    'ScH': 122,
    'SH': 123,
    'SiH': 124,
    'SiH2': 125,
    'SiH4': 126,
    'SiO': 127,
    'SiS': 129,
    'TiH': 130
}

atom_mass = {
    'H': 1.00794,
    'He': 4.002602,
    'Li': 6.941,
    'Be': 9.012182,
    'B': 10.811,
    'C': 12.0107,
    'N': 14.0067,
    'O': 15.9994,
    'F': 18.9984032,
    'Ne': 20.1797,
    'Na': 22.98976928,
    'Mg': 24.305,
    'Al': 26.9815386,
    'Si': 28.0855,
    'P': 30.973762,
    'S': 32.065,
    'Cl': 35.453,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
    'Sc': 44.955912,
    'Ti': 47.867,
    'V': 50.9415,
    'Cr': 51.9961,
    'Mn': 54.938045,
    'Fe': 55.845,
    'Co': 58.933195,
    'Ni': 58.6934,
    'Cu': 63.546,
    'Zn': 65.409,
    'Ga': 69.723,
    'Ge': 72.64,
    'As': 74.9216,
    'Se': 78.96,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.4678,
    'Sr': 87.62,
    'Y': 88.90585,
    'Zr': 91.224,
    'Nb': 92.90638,
    'Mo': 95.94,
    'Tc': 98.9063,
    'Ru': 101.07,
    'Rh': 102.9055,
    'Pd': 106.42,
    'Ag': 107.8682,
    'Cd': 112.411,
    'In': 114.818,
    'Sn': 118.71,
    'Sb': 121.76,
    'Te': 127.6,
    'I': 126.90447,
    'Xe': 131.293,
    'Cs': 132.9054519,
    'Ba': 137.327,
    'La': 138.90547,
    'Ce': 140.116,
    'Pr': 140.90465,
    'Nd': 144.242,
    'Pm': 146.9151,
    'Sm': 150.36,
    'Eu': 151.964,
    'Gd': 157.25,
    'Tb': 158.92535,
    'Dy': 162.5,
    'Ho': 164.93032,
    'Er': 167.259,
    'Tm': 168.93421,
    'Yb': 173.04,
    'Lu': 174.967,
    'Hf': 178.49,
    'Ta': 180.9479,
    'W': 183.84,
    'Re': 186.207,
    'Os': 190.23,
    'Ir': 192.217,
    'Pt': 195.084,
    'Au': 196.966569,
    'Hg': 200.59,
    'Tl': 204.3833,
    'Pb': 207.2,
    'Bi': 208.9804,
    'Po': 208.9824,
    'At': 209.9871,
    'Rn': 222.0176,
    'Fr': 223.0197,
    'Ra': 226.0254,
    'Ac': 227.0278,
    'Th': 232.03806,
    'Pa': 231.03588,
    'U': 238.02891,
    'Np': 237.0482,
    'Pu': 244.0642,
    'Am': 243.0614,
    'Cm': 247.0703,
    'Bk': 247.0703,
    'Cf': 251.0796,
    'Es': 252.0829,
    'Fm': 257.0951,
    'Md': 258.0951,
    'No': 259.1009,
    'Lr': 262,
    'Rf': 267,
    'Db': 268,
    'Sg': 271,
    'Bh': 270,
    'Hs': 269,
    'Mt': 278,
    'Ds': 281,
    'Rg': 281,
    'Cn': 285,
    'Nh': 284,
    'Fl': 289,
    'Mc': 289,
    'Lv': 292,
    'Ts': 294,
    'Og': 294,
    'ZERO': 0,
}

# Define the unit dictionary
unit = {
    'pc': 3.08567e16,        # m parsec
    'ly': 9.460730e15,       # m lightyear
    'AU': 1.49598e11,        # m astronomical unit
    'R_sun': 6.95700e8,      # m solar radius
    'R_jup': 7.1492e7,       # m nominal equatorial Jupiter radius (1 bar pressure level)
    'R_e': 6.371e6,          # m nominal Earth radius
    'd_H2': 2.827e-10,       # m molecular diameter of H2
    'M_sun': 1.989e30,       # kg solar mass
    'M_jup': 1.8982e27,      # kg Jupiter mass
    'M_e': 5.972e24,         # kg Earth mass
    'amu': 1.66054e-27,      # kg atomic mass unit
    'atm': 101325,           # Pa atmospheric pressure
}

# Define the const dictionary
pi = np.pi
const = {
    'k_B': 1.38065e-23,         # J K-1 Boltzmann constant
    'sig_B': 5.67037e-8,        # W m-2 K-4 Stephan Boltzmann constant
    'R': 8.31446,               # J mol-1 K-1 universal gas constant
    'G': 6.67430e-11,           # m3 kg-1 s-2 universal gravitational constant
    'eps_LJ': 59.7*5.67037e-8,  # J depth of the Lennard-Jones potential well for H2
    'c_p': 14300,               # J K-1 hydrogen specific heat
    'h': 6.62607,               # Js Planck's constant
    'hbar': 1.05457,            # Js
    'N_A': 6.02214e23,          # Avagadro's number
}

def mol_mass(M):
    """Calculate molecular mass in a.m.u. as long as the subscript of all
    elements in the chemical formula do not exceed 99. Input needs to be
    a string of case-correct chemical formula, e.g. CH4,
    """
    assert isinstance(M, str)
    assert M.isalnum()
    length = len(M)
    mass = 0
    while length > 0:
       if length == 1:
          # a single letter element
          mass = mass + atom_mass[M]
          break
       if length == 2:
          if M[1].isnumeric():
             mass = mass + atom_mass[M[0]]*int(M[1])
          elif M[1].upper() == M[1]:
             # two single letter elements
             mass += atom_mass[M[0]] + atom_mass[M[1]]
          else:
             mass +=  atom_mass[M]
          break
       else:
          if M[-1].isnumeric():
             # last char is number
             if M[-2].isnumeric():
                # last 2 char are numbers
                if M[-3].upper() == M[-3]:
                   # a single letter element
                   mass += atom_mass[M[-3]]*int(M[-2:])
                   M = M[:-3]
                else:
                   # a double letter element
                   mass += atom_mass[M[-4:-2]]*int(M[-2:])
                   M = M[:-4]
             else:
                # last char is number but second last char is string
                if M[-2].upper() == M[-2]:
                   # a single letter element
                   mass += atom_mass[M[-2]]*int(M[-1])
                   M = M[:-2]
                else:
                   # a double letter element
                   mass += atom_mass[M[-3:-1]]*int(M[-1])
                   M = M[:-3]
          else:
             # last char is string
             if M[-1].upper() == M[-1]:
                # a single letter element
                mass += atom_mass[M[-1]]
                M = M[:-1]
             else:
                # a double letter element
                mass += atom_mass[M[-2:]]
                M = M[:-2]
       length = len(M)
    return mass

def VERINT(X, Y, XIN):
    """
    C_DESCR:  two subroutines to perform vertical interpolation and integration
    C         of atmospheric profiles. They are stored together because the
    C        integration must implicitly assume the same interpolation scheme
    """

    if len(X)!=len(Y):
        raise('Interpolating from arrays of different sizes')
    if XIN>max(X) or XIN<min(X):
        print('Interpolating outside domain')

    if X[0]<X[-1]:
        for i in range(len(X)):
            if X[i]>XIN:
                I = i
                break
            else:
                I = len(X)
    else:
        for i in range(len(X)):
            if X[i]<XIN:
                I = i
                break
            else:
                I = len(X)
    if I == 0:
        I = 1

    if(X[I-2]!=X[I-1]):
        YOUT=Y[I-2]+(Y[I-1]-Y[I-2])*(XIN-X[I-2])/(X[I-1]-X[I-2])
    else:
        YOUT=Y[I-1]

    return YOUT
