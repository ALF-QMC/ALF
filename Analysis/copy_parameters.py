#!/usr/bin/env python3
"""Copies parameters from parameters file to hdf5 file.

"""
import sys
import os

import h5py
import f90nml

progdir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
    'Prog')
sys.path.append(progdir)
from parse_ham import parse


if len(sys.argv) > 1:
    filenames = sys.argv[1:]
else:
    filenames = ['./data.h5']

for filename in filenames:
    directory = os.path.dirname(filename)
    nml = f90nml.read(os.path.join(directory, 'parameters'))
    ham_name = nml['var_ham_name']['ham_name']

    default_parameters = parse(
        os.path.join(progdir, 'Hamiltonians',
                     'Hamiltonian_{}_smod.F90'.format(ham_name)))

    with h5py.File(filename, 'r+') as f:
        for nlist_name, nlist in default_parameters.items():
            nlist_name = nlist_name.lower()
            groupname = "parameters/{}".format(nlist_name)
            f.create_group(groupname)
            for par_name, par in nlist.items():
                par_name = par_name.lower()
                try:
                    val = nml[nlist_name][par_name]
                except KeyError:
                    val = par['value']

                if isinstance(par['value'], bool):
                    val = int(par['value'])
                f[groupname].attrs.create(par_name, val)
