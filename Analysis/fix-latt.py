#!/usr/bin/env python3
"""Fix a mix-up between `Norb` and `N_coord` in lattice."""
# pylint: disable=invalid-name

__author__ = "Jonas Schwab"
__copyright__ = "Copyright 2024, The ALF Project"
__license__ = "GPLv3"

from argparse import ArgumentParser
from glob import glob

import h5py


def _get_arg_parser():
    parser = ArgumentParser(
        description='Fix a mix-up between `Norb` and `N_coord` in lattice.',
        )
    parser.add_argument(
        'files', nargs='*',
        help='HDF5 files to fix. If empty, finds and fixes all '
            'files "data.h5", starting from the current working directory.')
    return parser

def _fix_latt(latt_grp):
    orbital_positions = []
    i = 0
    while True:
        i += 1
        try:
            orbital_positions.append(latt_grp.attrs[f"Orbital{i}"])
        except KeyError:
            break
    Norb = len(orbital_positions)

    if Norb != latt_grp.attrs["Norb"]:
        print("Fixing mixup between Norb and N_coord")
        print(f"Wrong values Norb={latt_grp.attrs['Norb']}, "
                f"N_coord={latt_grp.attrs['N_coord']}")
        print(f"Correct values Norb={latt_grp.attrs['N_coord']}, "
                f"N_coord={latt_grp.attrs['Norb']}")
        if Norb != latt_grp.attrs['N_coord']:
            raise RuntimeError("Non-recoverable mixup with Norb/N_coord")
        latt_grp.attrs["N_coord"] = latt_grp.attrs["Norb"]
        latt_grp.attrs["Norb"] = Norb


def _fix_loc(hdf5_loc):
    for i in hdf5_loc:
        if isinstance(i, str):
            print(i)
            if i == 'lattice':
                _fix_latt(hdf5_loc[i])
            else:
                _fix_loc(hdf5_loc[i])


if __name__ == '__main__':
    args = _get_arg_parser().parse_args()

    files = args.files if args.files else glob("**/data.h5", recursive=True)
    print('foo', files)
    for fname in files:
        print(f'Fixing "{fname}"')
        with h5py.File(fname, 'r+') as f:
            _fix_loc(f)
