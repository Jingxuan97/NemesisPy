#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json
"""
Routine to turn a RADTRAN style gas data file into a Python dictiondary
of the form:
gas_info = {
    "1": {
        "name": "H2O",
        "isotope": {
            "1": {
                "abun": 0.997317,
                "mass": 18.0106,
                "id": 161,
                "partition": [
                    -3.46229,
                    0.251974,
                    0.00143133,
                    -8.45955e-07
                ]
            },
            "2" :{...}
        },
        "mmw": 18.01529586265
    },
    "2": {...},
    "3": {...},
}

Note: need to remove the comments in the begining of the file and make sure
that there are no spaces in the isotop id definition, i.e.
dont want (  1) or ( 23). This is already done in 'gasinforef.txt'.
"""
gas_data_file = 'gasinforef.txt'
gas_info = {}
numer_of_molecules = 130
def rmsp(s): # remove leading blank spaces
    while s[0] == ' ':
        s = s[1:]
    return s

with open(gas_data_file) as fp:
    for i in range(numer_of_molecules):
        line = fp.readline()

        if line[1] == '-' or line[2] == '-':
            line =  rmsp(fp.readline())
            # print('new molecule')

            rank = line.split()[0]
            # print('molecule rank',int(rank))
            line =  rmsp(fp.readline())

            name = line.split()[0]
            # print('name', str(name))

            line =  rmsp(fp.readline())
            n_iso = line.split()[0]
            # print('n_iso', int(n_iso))
            n_iso = int(n_iso)

            mol = {'name': '{}'.format(name), 'isotope':{}}
            # print(mol)

            for i in range(n_iso):
                line =  rmsp(fp.readline())
                # print('isotope', i+1)
                isorank, abun, mass, isoid = line.split()[0], line.split()[1], line.split()[2], line.split()[3]
                isoid = isoid[1:-1]
                # print('isorank, abun, mass, isoid',isorank, abun, mass, isoid)

                line =  rmsp(fp.readline())
                partition = line.split()
                for i in range(len(partition)):
                    partition[i] = float(partition[i])
                # print('partition', partition)

                mol['isotope']['{}'.format(isorank)] = {'abun': float(abun),
                                                        'mass': float(mass),
                                                        'id': int(isoid),
                                                        'partition': partition}
            mmw = 0
            for i in range(n_iso):
                mmw += mol['isotope']['{}'.format(i+1)]['mass'] \
                    * mol['isotope']['{}'.format(i+1)]['abun']
            mol['mmw'] = mmw

            rank = int(rank)
            # print(type(rank))
            gas_info[(int(rank))] = mol


print('------------------------------------------------')
print('------------------------------------------------')

# print(json.dumps(gas_info, indent=4, sort_keys=False))
f = open('gas_info.py', 'w')
f.write('gas_info = ')
f.write(json.dumps(gas_info, indent=4, sort_keys=False))
f.close()

